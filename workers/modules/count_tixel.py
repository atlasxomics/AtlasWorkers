# --------------------------------------- #
# Small library of functions utilized for segmenting
# and counting stained Nuclei from a tissue image.
#
# Libraries needed for install: Opencv, numpy, skimage, cellpose scikit-image
#
# company: AtlasXomics
# email: jonahs@atlasxomics.com
# ---------------------------------------- #
from cellpose import models
from cellpose.io import imread
import numpy as np
import random
import cv2
from skimage.feature import peak_local_max
import csv


class NucleiCounter():
    cellpose_model = models.Cellpose(gpu = False, model_type="nuclei")

    def __init__(self,
        min_area_in_tixel = 50,
        nuclei_size = 280,
        nuclei_size_sf = 1.6,
        threshold = 70,
        blur_kernel_size = 5,
        number_blur_iters = 3,
        min_distance_peak_local_max = 4,
        opening_kernel_size = 3,
        opening_kernel_iterations = 1,
        dialation_iterations = 3,
        channels = [0, 0],
        model = cellpose_model
    ):
        self.min_area_in_tixel =  min_area_in_tixel
        self.nuclei_size = nuclei_size
        self.nuclei_size_thresh = nuclei_size_sf * self.nuclei_size
        self.threshold = threshold
        self.blur_kernel_size = blur_kernel_size
        self.number_blur_iters = number_blur_iters
        self.min_distance_peak_local_max = min_distance_peak_local_max
        self.opening_kernel_size = opening_kernel_size
        self.opening_kernel_iterations = opening_kernel_iterations
        self.dialation_iterations = dialation_iterations
        self.channels = channels
        self.cellpose_model = model
        
    def get_tixel_rand(self, img, width):
        """ Returns a tuple of the new subset of the larger image 
            and the x and y coordinates used to create it
            
            img: the full image being worked with
            width: The width of each tixel in pixels
        """
        "Returns a random section of the provided img that is 2x provided width in height and width"
        w = img.shape[1]
        h = img.shape[0]
        x = int(round((random.random() * (w - (2 * width))) + width))
        y = int(round((random.random() * (h - (2 * width))) + width))
        tix = img[y - width: y + width, x - width: x + width].copy()
        return (tix, x, y)


    def get_tixel_coord(self, img, x, y, width):
        """Returns the subsetted image based on given coordinates and tixel width
            # img: the full image being worked with
            x: x-coord
            y: y-coord
            width: the width of each tixel in pixels
        """
        t = img[y - width: y + width, x - width: x + width].copy()
        return t

    #
    def count_tixel(self, tixel, width, path, inx, sizes, channel, methods, diameter = None):
        """Takes in the tixel image and the known width of the tixel and makes predictions about number of nuclei within it.
        Three methods possible: 0: Watershed 1: Size segmentation 2: Cellpose Nuclei
        tixel: The subsetted image that is 2x tixel width and height at the desired location on image
        width: Tixel width in pixels
        path: Where to store the predicted image
        inx: The index of this tixel prediction. Used to keep track of count
        sizes: List of area of tixels predicted
        channel: The channel of the image the color of nuclei are on
        method: 0: Watershed 1: Threshold Area Segmentation 2: Cellpose Prediction
        """

        image_method_preds = {}

        blue = tixel[:,:, channel]
        t_orig = tixel.copy()
        t = blue.copy()

        if 2 in methods:
            # for b in range(self.number_blur_iters):
            #     t = cv2.blur(t, (self.blur_kernel_size, self.blur_kernel_size))
            pred = self.pred_cellpose_segmentation(t, tixel, diameter, inx, width, path)
            image_method_preds["cellpose"] = pred 

        if 0 in methods or 1 in methods:
            for b in range(self.number_blur_iters):
                t = cv2.blur(t, (self.blur_kernel_size, self.blur_kernel_size))
                t_orig = cv2.blur(t_orig, (self.blur_kernel_size, self.blur_kernel_size))

            ret, thresh = cv2.threshold(t, self.threshold, 255, cv2.THRESH_BINARY)
            kernel = np.ones((self.opening_kernel_size,self.opening_kernel_size), np.uint8)
            opening = cv2.morphologyEx(thresh, cv2.MORPH_OPEN, kernel, iterations=self.opening_kernel_iterations)
            if 0 in methods:
                pred = self.pred_watershed(t_orig, tixel, opening, kernel, path, inx, width, sizes)
                image_method_preds["watershed"] = pred
            if 1 in methods:
                pred = self.pred_threshold_segmentation(tixel, thresh, width, path, inx)
                image_method_preds["thresh_seg"] = pred

        return image_method_preds



    def pred_watershed(self, t_orig, tixel_nonblur, opening,kernel, path, inx, width, sizes):
        """Prediction of nuclei within the given image using an implementation of the watershed algorithm as well as size segmentation
            t_orig: The 2x width and height image around the coordinate of interest which has been blurred
            tixel_nonblur: The same image as t_orig but it has not been blurred
            opening: Same image after having opening (erosion followed by dialation) preformed
            path: where to store images
            inx: Index of prediction
            width: Width of tixel in pixels
            sizes: List of predicted nuclei areas
        
        """
        "Returns the number of predicted nuclei"

        "Setup and determination of watershed image"
        surebg = cv2.dilate(opening, kernel, iterations=self.dialation_iterations)
        D = cv2.distanceTransform(opening, cv2.DIST_L2, 5)
        maxi = peak_local_max(D, min_distance=self.min_distance_peak_local_max, exclude_border=False)
        surefg = np.zeros(D.shape, dtype=np.uint8)
        surefg[tuple(maxi.T)] = 255
        unknown = cv2.subtract(surebg, surefg)
        ret, markers = cv2.connectedComponents(surefg)
        markers = markers + 1
        markers[unknown == 255] = 0
        water = cv2.watershed(t_orig, markers)

        "Final prediction taking into account nuclei size"
        (coord_1, coord_2) = self.get_inner_width(t_orig, width)
        tixel_water = water[coord_1: coord_2, coord_1: coord_2]
        vals_total = np.unique(water)
        vals = np.unique(tixel_water)
        counted_colors = []
        pred = 0
        total_area = 0
        for k in range(len(vals_total)):
            val = vals_total[k]
            if val != 1 and val != -1:
                size = np.count_nonzero(water == val)
                sizes.append(size)
                if val in vals:
                    size_tixel = np.count_nonzero(tixel_water == val)
                    total_area += size_tixel
                    if size_tixel >  self.min_area_in_tixel:
                        guess = 1
                        counted_colors.append(val)
                        if size_tixel > self.nuclei_size:
                            guess = int(size_tixel / self.nuclei_size)
                        pred += guess

        if path is not None:
            self.save_prediction_image(water,counted_colors, tixel_nonblur, width, inx, path, name_method="watershed")

        return pred

    def save_prediction_image(self, highlighted_image, color_values, original_image, tixel_width, inx, path, name_method = "highlighted"):
        """Takes in the `image` with markings of segementation and saves the result in specified spot.
            Saves a version that highlights the different blobs and one of the original
        """
        original_image = cv2.cvtColor(original_image, cv2.COLOR_BGR2RGB)
        original_image_copy = original_image.copy()
        original_image_copy2 = original_image.copy()
        (coord_1, coord_2) = self.get_inner_width(original_image_copy, tixel_width)
        for color_val in color_values:
            c1 = (random.random() * 255)
            c2 = (random.random() * 255)
            c3 = (random.random() * 255)

            original_image_copy[highlighted_image == color_val] = [c1, c2, c3]


        cv2.rectangle(original_image_copy, (coord_1, coord_1), (coord_2, coord_2), (0, 255, 255), 1)
        cv2.rectangle(original_image_copy2, (coord_1, coord_1), (coord_2, coord_2), (0, 255, 255), 1)


        path_original = f"{path}/tixel_{str(inx)}.tif"
        path_highlighted = f"{path}/tixel_{name_method}_{str(inx)}.tif"
        cv2.imwrite(path_original,original_image_copy2)
        cv2.imwrite(path_highlighted, original_image_copy)



    def pred_cellpose_segmentation(self, tixel, tixel_original, diameter, inx, tixel_width, path = None):
        """Uses the Cellpose nuclei segmentation model to predict the number of nuclei in a tixel of an image"""
        masks, flows, styles, diams = self.cellpose_model.eval(tixel, diameter=diameter, channels=self.channels)
        (coord1, coord2) = self.get_inner_width(tixel, tixel_width)
        tixel_only_mask = masks[coord1: coord2, coord1: coord2]
        total_colors = np.unique(masks)
        colors_count = []
        predicted_colors = [color for color in total_colors if color > 0]
        count = 0
        for color in predicted_colors:
            tixel_area = np.count_nonzero(tixel_only_mask == color)
            if tixel_area >  self.min_area_in_tixel:
                colors_count.append(color)
                count += 1
        if path:
            self.save_prediction_image(masks, colors_count, tixel_original, tixel_width, inx, path, name_method = "cellpose")
        return count

    def pred_threshold_segmentation(self, t_orig, opening, width, path, inx):
        "Takes a binary threshold of the image and predicts count based on size"
        (coord_1, coord_2) = self.get_inner_width(t_orig, width)
        open_disp = opening.copy()
        disp_orig = t_orig.copy()
        inner_tixel = opening[coord_1: coord_2, coord_1: coord_2]

        in_area = (np.sum(np.where(inner_tixel == 255, 1, 0)))
        count = int(round(in_area / self.nuclei_size_thresh))
        if path is not None:
            open_disp = cv2.rectangle(open_disp, (coord_1, coord_1), (coord_2, coord_2), 125, 1)
            disp_orig = cv2.rectangle(disp_orig, (coord_1, coord_1), (coord_2, coord_2), (0, 255, 0), 1)
            disp_orig = cv2.cvtColor(disp_orig, cv2.COLOR_BGR2RGB)
            cv2.imwrite(f"{path}/bw_{inx}.tif", open_disp)
            cv2.imwrite(f"{path}/orig_{inx}.tif", disp_orig)

        return count

    def pred_thresh_segmentation_whole_image(self, thresholded):
        "Takes a binary threshold of the image and predicts count based on size"
        # for b in range(self.number_blur_iters):
            # thresholded = cv2.blur(thresholded, (self.blur_kernel_size, self.blur_kernel_size))
            #  = cv2.blur(t_orig, (self.blur_kernel_size, self.blur_kernel_size))

        # kernel = np.ones((self.opening_kernel_size, self.opening_kernel_size), np.uint8)
        # opening = cv2.morphologyEx(thresholded, cv2.MORPH_OPEN, kernel, iterations=self.opening_kernel_iterations)
        count = (np.sum(np.where(thresholded == 255, 1, 0)) // self.nuclei_size_thresh)
        return count
        
    
    def get_inner_width(self, img, width):
        """Returns the coordinates defining the top left and bottom right of a inner tixel based on tixel width"""
        coord = img.shape[0]
        c1 = (coord // 2) - (width // 2)
        c2 = (coord // 2) + (width // 2)
        return (c1, c2)

    def in_contours(self, point, contours, in_contours = True):
        "Method to test whether a randomly generated point is within any of the listed contours."
        for contour in contours:
            if cv2.pointPolygonTest(contour, point, False) >= 0:
                return in_contours
        return (not in_contours )

    def survey_image(self, img, width_pixels, channel = 0, path = None, n = 100,  predictions = {}, tixel_sizes = [], contours = None, methods = [0], diameter = None, in_contours = True):
        "Method to preform an image study. Generating random subsets of the larger image and predicting counts on a tixel within the subset."
        i = 0
        while i < n:
            (t, x, y) = self.generate_sub_image(img, width_pixels, contours=contours, in_contours=in_contours)
            pred = self.count_tixel(t, width_pixels, path, i, tixel_sizes, channel, methods = methods, diameter = diameter)

            predictions[i] = {"x": x, "y": y, "counts": pred}
            i += 1
        return predictions
    
    def generate_sub_image(self, img, width_pixels, contours = None, in_contours = True):
        """Generates a random subset image with a midpoint that falls within the defined space.

        Args:
            img (np array): The overall image to sample from
            width_pixels (_type_): Width in pixels of the hypothetical tixel
            contours (list[numpy arrays], optional): List of contour shapes used to define proper space. Defaults to None.
            in_contours (bool, optional): Whether to used points within the contours (True) or outside (False). Defaults to True.

        Returns:
            numpy array: Sub image that is in the user defined space
        """
        while True:
            (t, x, y) = self.get_tixel_rand(img, width_pixels)
            if contours is not None:
                correct_place = self.in_contours((x, y), contours=contours, in_contours=in_contours)
                if correct_place:
                    break
            else:
                break
        return (t, x, y)
                    
    def find_nuclei_sizes(self, image_grayscale, sizes = []):
        """Returns the size of each nuclei identified in the image.

        Args:
            image_grayscale (np.array): gray scale image being used to count on
            sizes (list, optional): List where sizes are appended. Defaults to [].

        Returns:
            list: Size of each nuclei, in form {color: area}
        """
        masks, _, _, diameter = self.cellpose_model.eval(image_grayscale, diameter = None, channels=self.channels)
        total_colors = [color for color in np.unique(masks) if color > 0]
        sizes = []
        for color in total_colors:
            size = np.count_nonzero(masks == color)
            res = {color: size}
            sizes.append(res)
        return sizes, diameter
    
    def load_contours(self, filename):
        """ Used to load countours defined in a csv file in as numpy arrays.
        

        Args:
            filename (string): path to the file containing the contours
            
        Return: list of numpy arrays
        """
        def end_contour(current, lis):
            np_curr = np.array(current)
            lis.append(np_curr)
        contents = open(filename, "r")
        csv_reader = csv.reader(contents, delimiter = ",")
        contour_lis = []
        current_contour = []
        for inx, row in enumerate(csv_reader):
            if inx == 0: continue
            if row[0][0] == "C":
                if current_contour:
                    end_contour(current_contour, contour_lis)
                    current_contour = []
                continue
            x = int(row[0])
            y = int(row[1])
            inner = np.array([x, y])
            outer = np.array([inner])
            current_contour.append(outer)
        
        end_contour(current_contour, contour_lis)
        return contour_lis
    
    def get_mask_image_cellpose(self, img, roi_boundary1,roi_boundary2, diameter = None, channels = [0,0]):
        """Returns a mask image of the region of interest.

        Args:
            img (np.array): The image to mask
            roi_boundary1 (tuple): The top left corner of the region of interest
            roi_bounday2 (tuple): The bottom right corner of the region of interest
        """
        
        x1, y1 = roi_boundary1
        x2, y2 = roi_boundary2
        
        sub_img =img[y1:y2, x1:x2]
        masks, _, _, diameter = self.cellpose_model.eval(sub_img, diameter = None, channels=channels)
        return masks
        
        
        
        
    def generate_tixel_mapping_cellpose(self, img, tixel_width_pixels, tissue_pos_file):
        """
        Method to append the counts of nuclei in a tixel to the tissue position file.
        
        img np.array: Registered image to count on.
        tixel_width_pixels int: The width of the tixel in pixels.
        tissue_pos_file Pandas Dataframe: The loaded tissue position file. Contains barcode: string, on_tissue: int ,row_inx: int, col_inx: int, y coord: int x coord: int.
        
        return: Pandas Dataframe: The tissue position file with the cell count appended.
        """
        mask, _, _, diam = self.cellpose_model.eval(img, diameter = None, channels = [0,0])
        
        tissue_pos_file_cellcount = tissue_pos_file.copy()
        tissue_pos_file_cellcount["cell_count"] = 0
        for index, row in tissue_pos_file.iterrows():
            
            x_coord = row["x_coord"]
            y_coord = row["y_coord"]
            
            sub_tixel = self.get_tixel_coord(mask, x_coord, y_coord, (tixel_width_pixels // 2))
            # plt.imshow(sub_tixel)
            # plt.show()
            count = self.count_cells_in_tixel(sub_tixel)
            
            tissue_pos_file_cellcount.at[index, "cell_count"] = count
        
        # write_file.close()
        return tissue_pos_file_cellcount, count
            
        
    def count_cells_in_tixel(self, tixel_mask):
        """Counts the number of cells in a tixel.

        Args:
            tixel (np.array): The tixel to count cells in. This should be a masked image where each unique pixel intensity > 0 represents a cell.

        Returns:
            int: The number of cells in the tixel
        """
        total_colors = [color for color in np.unique(tixel_mask) if color > 0]
        counts = 0
        for color in total_colors:
            size = np.count_nonzero(tixel_mask == color)
            if size > self.min_area_in_tixel:
                counts += 1
        return counts