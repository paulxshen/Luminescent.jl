# # -*- coding: utf-8 -*-
# """Convert an image file to a GDS file
# """

# import cv2
# import numpy as np
# import gdspy

# import argparse


# def minmax(v):
#     if v > 255:
#         v = 255
#     if v < 0:
#         v = 0
#     return v


# def main(fileName, sizeOfTheCell, layerNum, isDither, scale):
#     """Convert an image file (fileName) to a GDS file
#     """
#     print("Converting an image file to a GDS file..")
#     # Read an image file
#     img = cv2.resize(cv2.imread(fileName), dsize=None, fx=scale, fy=scale)

#     width = img.shape[1]
#     height = img.shape[0]
#     print("width:{0}".format(width))
#     print("height:{0}".format(height))

#     # Convert an image to grayscale one
#     gray = cv2.cvtColor(img, cv2.COLOR_RGB2GRAY)

#     if (isDither):
#         # Floyd–Steinberg dithering
#         # https://en.wikipedia.org/wiki/Floyd%E2%80%93Steinberg_dithering
#         for y in range(0, height-1):
#             for x in range(1, width-1):
#                 old_p = gray[y, x]
#                 new_p = np.round(old_p/255.0) * 255
#                 gray[y, x] = new_p
#                 error_p = old_p - new_p
#                 gray[y, x+1] = minmax(gray[y, x+1] + error_p * 7 / 16.0)
#                 gray[y+1, x-1] = minmax(gray[y+1, x-1] + error_p * 3 / 16.0)
#                 gray[y+1, x] = minmax(gray[y+1, x] + error_p * 5 / 16.0)
#                 gray[y+1, x+1] = minmax(gray[y+1, x+1] + error_p * 1 / 16.0)

#         ret, binaryImage = cv2.threshold(gray, 0, 255, cv2.THRESH_OTSU)
#     else:
#         ret, binaryImage = cv2.threshold(gray, 0, 255, cv2.THRESH_OTSU)

#     # Fill orthological corner
#     for x in range(width - 1):
#         for y in range(height - 1):
#             if binaryImage.item(y, x) == 0 and binaryImage.item(y + 1, x) == 255 \
#                     and binaryImage.item(y, x + 1) == 255 and binaryImage.item(y + 1, x + 1) == 0:
#                 # binaryImage.itemset((y + 1, x), 0)
#                 binaryImage[y + 1, x] = 0
#             elif binaryImage.item(y, x) == 255 and binaryImage.item(y + 1, x) == 0 \
#                     and binaryImage.item(y, x + 1) == 0 and binaryImage.item(y + 1, x + 1) == 255:
#                 # binaryImage.itemset((y + 1, x + 1), 0)
#                 binaryImage[y + 1, x+1] = 0

#     # Output image.bmp
#     cv2.imwrite("image.bmp", binaryImage)

#     # The GDSII file is called a library, which contains multiple cells.
#     lib = gdspy.GdsLibrary()
#     gdspy.current_library = gdspy.GdsLibrary()

#     # Geometry must be placed in cells.
#     unitCell = lib.new_cell('CELL')
#     square = gdspy.Rectangle((0.0, 0.0), (1.0, 1.0), layer=(int)(layerNum))
#     unitCell.add(square)

#     grid = lib.new_cell("GRID")

#     for x in range(width):
#         for y in range(height):
#             if binaryImage.item(y, x) == 0:
#                 # print("({0}, {1}) is black".format(x, y))
#                 cell = gdspy.CellReference(
#                     unitCell, origin=(x, height - y - 1))
#                 grid.add(cell)

#     scaledGrid = gdspy.CellReference(
#         grid, origin=(0, 0), magnification=(float)(sizeOfTheCell))

#     # Add the top-cell to a layout and save
#     top = lib.new_cell("TOP")
#     top.add(scaledGrid)
#     s = fileName.split(".")
#     s = (".").join(s[:-1])
#     lib.write_gds(f"{s}.gds")


# if __name__ == "__main__":
#     parser = argparse.ArgumentParser()
#     parser.add_argument('fileName', type=str,
#                         help='name of the input image file')
#     parser.add_argument('sizeOfTheCell', type=float,
#                         help='size of the unit-cells (minimum width and space) [um]')
#     parser.add_argument('layerNum', type=int,
#                         help='layer number of the output GDSII file')
#     parser.add_argument('--scale', default=1.0, type=float, help='scale')
#     parser.add_argument('-d', action='store_true',
#                         help='Floyd–Steinberg dithering')
#     args = parser.parse_args()

#     main(args.fileName, args.sizeOfTheCell, args.layerNum, args.d, args.scale)
