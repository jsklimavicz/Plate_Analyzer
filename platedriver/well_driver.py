# platedriver/well_driver.py
# This file is a component of Plate_Analyzer, which can count Drosophila
# L1 larvae, classify them as alive or dead, and determine dose-repsonse
# information based on live/dead count data. 

# Copyright (C) 2021 James Klimavicz

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.


import well_prep as wpp
import plate_reflected as Plate


# folder = '/home/jklimavicz/Documents/Merlin_Images/210526_132949_05262021 plate 1/210526_132949_Plate 1/'

folders = []
# folders.append('/home/jklimavicz/Documents/Merlin_Images/210611_130231_06112021_Roadmap/210611_130231_Plate 1')
# folders.append('/home/jklimavicz/Documents/Merlin_Images/210610_133258_06102021_WaQuita mal test/210610_133258_Plate 1')
# folders.append('/home/jklimavicz/Documents/Merlin_Images/210528_115911_05282021_colors/210528_115911_Plate 1')
# folders.append('/home/jklimavicz/Documents/Merlin_Images/210526_132949_05262021 plate 1/210526_132949_Plate 1')
# folders.append('/home/jklimavicz/Documents/Merlin_Images/210528_152349_05282021_Roadmap1/210528_152349_Plate 1')
# folders.append('/home/jklimavicz/Documents/Merlin_Images/210611_130231_06112021_Roadmap/210611_130231_Plate 1')
folders.append('/home/jklimavicz/Documents/Merlin_Images/210623_143228_06232021_Roadmap 3/210623_143228_Plate 1')
folders.append('/home/jklimavicz/Documents/Merlin_Images/210624_141637_06242021_Roadmap/210624_141637_Plate 1')



save_folder = '/home/jklimavicz/Documents/Merlin_Images/Merlin_Annotation/new/train_orig'
mask_folder = '/home/jklimavicz/Documents/Merlin_Images/Merlin_Annotation/new/mask'

# folders = ['/home/jklimavicz/Documents/Merlin_Images/210528_152349_05282021_Roadmap1/210528_152349_Plate 1']

for folder in folders:
	print(folder)
	plate = Plate.Plate(folder)
	plate.run_images(init_only=False)
	plate.save_images(save_folder, mask_folder=mask_folder)


# fig, ax = plt.subplots(figsize=(12, 9))
# ax.imshow(plate.full_plate)
# plt.show()
exit()



