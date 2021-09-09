Analysis of flourescently labeled cholesterol sensor images 
================

## Introduction
Analysis was performed in Python. Images were segmented using Cellpose Python library. To select high-quality images the cytoplasm marker channel was used with Laplace filtering. We used high-quality images as input of the Cellpose model with parameter channel set to greyscale and cell diameter greater than 200 pixels. After identifying cell boundaries, we applied binary erosion with default structure and 10 iterations to determine cytoplasm boundary, or binary dilation with default structure and 5 iterations to determine PM outer boundary. The boundary of PM was determined by subtracting the cytoplasm boundary from the outer boundary. We calculated the log2 ratio of the mean PM and mean intracellular D4H fluorescence intensities for each cell in the D4H channel to examine the changes of plasma membrane cholesterol distribution. 

### Creating metadata for images and filtering for the analysis

``` py
# Import required packages
import numpy as np
import pandas as pd
import os
import sys

from scipy.stats import spearmanr
from scipy.stats import pearsonr

from skimage.filters import laplace
from cellpose import utils, io, models, transforms

import importlib
import importlib.util

# Import experiment related constant values (path, meta filename, etc.) as parameter in command line
spec = importlib.util.spec_from_file_location("consts", "./consts/"+sys.argv[1])
c = importlib.util.module_from_spec(spec)
spec.loader.exec_module(c)
print('Consts read in from', sys.argv[1])

# Read in data
timepoints = os.listdir(c.raw)
inames = []
for tp in timepoints:
    inames += [c.raw + '/' + tp + '/' + x for x in os.listdir('%s/%s' % (c.raw, tp))
               if (('thumb' not in x) & ('.tif' in x))]



meta = pd.DataFrame(index=range(len(inames)))
meta['Raw'] = inames
meta['TimePoint'] = meta['Raw'].apply(
    lambda x: int(x.split('/')[c.pos+1].split('_')[1]))
meta['Well'] = meta['Raw'].apply(lambda x: x.split('/')[c.pos+2].split('_')[1])
meta['Region'] = meta['Raw'].apply(lambda x: x.split('/')[c.pos+2].split('_')[2])
meta['WaveLength'] = meta['Raw'].apply(
    lambda x: x.split('/')[c.pos+2].split('_')[3][:2])
meta['Laplace'] = 0.0
meta['Laplace_norm'] = 0.0
meta['Laplace_std'] = 0.0
meta['Laplace_std_norm'] = 0.0
meta = meta.sort_values(['TimePoint', 'Well', 'Region', 'WaveLength'])
meta.index = range(len(meta.index))

meta['Mean'] = 0.0
meta['Median'] = 0.0
meta['Std'] = 0.0
meta['Min'] = 0.0
meta['Max'] = 0.0
meta = meta[['TimePoint', 'Well', 'Region', 'WaveLength',
             'Mean', 'Median', 'Std', 'Max', 'Min', 'Laplace', 'Laplace_norm', 'Laplace_std', 'Laplace_std_norm',
             'Raw']]

for i in meta.index:
    if i % 1000 == 0:
        print('Done for %i images!' % i)
    path_old = meta.loc[i, 'Raw']
    img = io.imread(path_old)
    img_n = transforms.normalize99(img)
    lplc = laplace(img)
    lplc_norm = laplace(img_n)
    meta.loc[i, ['Mean', 'Median', 'Std', 'Max', 'Min',
                 'Laplace', 'Laplace_norm', 'Laplace_std', 'Laplace_std_norm']] = (img.mean(), np.median(img), img.std(), img.max(), img.min(),
                                                lplc.max(), lplc_norm.max(), lplc.std(), lplc_norm.std())

# Save metadata to expreiments
meta.to_csv(f'{c.path}/{c.meta_name}.csv')

# Filter metadata
# Filter well (plasmamembrane or cholesteol sensor) and  laplace_norm values
fil = (meta['WaveLength'] == 'w2') & (
    meta['Well'].apply(lambda x: x[0] in c.wells_to_filter))
meta = meta[fil]
fil = meta[c.filter_column] < c.laplace_th
meta = meta[fil]

meta.to_csv(f'{c.path}/{c.meta_filtered_name}.csv')
print('Filtered metadata saved')
```

Imported constant files contain path of raw images,metadata and filtered metadata, masks and results filename, threshold value for laplace filtering, column to filter, wells to filter. 
``` py

# Example const.py file for an experimnet
raw = '../../../ubuntu/images/SARS-CoV-2/20210128a/2021-01-28/1410/'
pos = raw.count('/')
path = '../exp2'

laplace_th = 0.15
th = str(laplace_th).replace('.', '')
wells_to_filter = ['C', 'D']
filter_column = 'Laplace_std_norm'

meta_name = 'exp2a_images_meta'
meta_filtered_name = f'exp2a_meta_laplace_std{th}_w2'
masks_folder = f'{path}/masks_std{th}_a'
results_name = f'exp2a_results_lp_std{th}_mp'
```



### Calculating masks using Cellpose library

``` py
# Import required packages
import numpy as np
import pandas as pd

import os
import os.path
import sys

from scipy.ndimage.morphology import binary_dilation, binary_erosion

from cellpose import utils, io, models, transforms

import importlib
import importlib.util

# Import experiment related constant values (path, meta filename, etc.) as parameter in command line
spec = importlib.util.spec_from_file_location("consts", "./consts/"+sys.argv[1])
c = importlib.util.module_from_spec(spec)
spec.loader.exec_module(c)
print('Consts read in from', sys.argv[1])

# Fit Cellpose model
print('Fit model')
model = models.Cellpose(gpu=False, model_type='cyto')
channels = [0, 0]

# Read in filtered metadata
# filtered with Laplace threshold and wavelength
meta_filtered = pd.read_csv(f'{c.path}/{c.meta_filtered_name}.csv', index_col = 0)

# Fit cellpose model
print('Fit model')
model = models.Cellpose(gpu=False, model_type='cyto')
channels = [0, 0]
print('')

# Calulate masks for every image and save them
print('Calculate masks')
for i in meta_filtered.index:
    if (os.path.isfile(f'{c.masks_folder}/{i}_cp_masks.png')):
        print(f'{c.masks_folder}/{i}_cp_masks.png exists')
        continue
    else:
        img = io.imread(meta_filtered.loc[i, 'Raw'])
        masks, flows, styles, diams = model.eval(
            img, diameter=200, channels=channels, net_avg=False)
        io.save_to_png(img, masks, flows, f'{c.masks_folder}/{i}.png')
print('Done')
```

### Calculating cytoplasm and membrane intensities (PM/IC ratio) of cells 

``` py
# import required packages
import numpy as np
import pandas as pd
import os
import sys

from scipy.ndimage.morphology import binary_dilation, binary_erosion
from skimage.filters import laplace

from cellpose import utils, io, models, transforms

from multiprocessing import Pool
import time

import importlib
import importlib.util

# Import experiment related constant values (path, meta filename, etc.) as parameter in command line
spec = importlib.util.spec_from_file_location("consts", "./consts/"+sys.argv[1])
c = importlib.util.module_from_spec(spec)
spec.loader.exec_module(c)
print('Consts read in from', sys.argv[1])

spec = importlib.util.spec_from_file_location("mapping", "./consts/mapping.py")
mapping = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mapping)
print('Wellmap read in from mapping.py')

# Functions

# Calculating cell outer and inner boundaries of cell using binary erosion and dilatation
def extend_mask(m):
    m1 = m.copy()
    m2 = m.copy()
    for i in range(5):
        m1 = binary_dilation(m1)
    for i in range(10):
        m2 = binary_erosion(m2)
    return m1*1, m2*1

# Calculating cytoplasm and membrane intensities of all cells on an image 
def calculate_cyto_membrane_intensity(i):
    results = pd.DataFrame(
        columns=['ID_w1', 'ID_w2', 'Mask', 'Cyto', 'Membrane'])

    img2 = io.imread(meta.loc[i-1, 'Raw'])

    mask = io.imread(f'{c.masks_folder}/{i}_cp_masks.png')
    cyto_ = []
    memb_ = []
    for j in range(1, mask.max()+1):
        m = mask == j
        m1, m2 = extend_mask(m)
        cyto = m2 == 1
        memb = (m1-m2) == 1
        cyto = (cyto*img2).sum() / cyto.sum()
        memb = (memb*img2).sum() / memb.sum()
        cyto_.append(cyto)
        memb_.append(memb)
    temp = pd.DataFrame(index=range(len(cyto_)), columns=results.columns)
    temp['ID_w1'] = i-1
    temp['ID_w2'] = i
    temp['Mask'] = list(range(1, mask.max()+1))
    temp['Cyto'] = cyto_
    temp['Membrane'] = memb_
    results = pd.concat([results, temp])
    # results.index = range(len(results.index))
    return results

print('Read in masks and calculate results')

# Read in data

meta = pd.read_csv(f'{c.path}/{c.meta_filtered_name}.csv', header=0, index_col=0)
fnames = list(meta.index.sort_values())

meta = pd.read_csv(f'{c.path}/{c.meta_name}.csv', header=0, index_col=0)

print('Number of images:', len(fnames))

# Calculating the intensities of all images (using multiprocessing pool)
p = Pool()
start_time = time.time()
result_list = p.map(calculate_cyto_membrane_intensity, fnames)
p.close()
p.join()

# Printing the time of the process
end_time = time.time()-start_time
print(str(round(end_time/60))+' min')

# Calculating ratio of intensities and save results
results = pd.concat(result_list)
results = results.reset_index(drop=True)

file = f'{c.path}/{c.results_name}.csv'
print(f'Save results table to {file}')
results.to_csv(file)

complete_file = f'{c.path}/{c.results_name}_complete.csv'
print(f'Merge results and meta, calculate ratio and save it to {complete_file}')
results = pd.merge(results, meta, left_on='ID_w1', right_index=True)
results['Ratio'] = results['Membrane'] / results['Cyto']
results['Drug'] = results['Well'].map(mapping.wellmap)
results['logRatio'] = np.log2(results['Ratio'])
results['WellRegion'] = results['Well'] + results['Region']

# Save results
results.to_csv(complete_file)
print('Done')

```
