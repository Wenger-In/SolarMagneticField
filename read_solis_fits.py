from astropy.io import fits
import numpy as np

# 打开.fits文件
fits_file = fits.open('E:/Research/Data/SOLIS/kbv7g100506t0941c2096_000_int-mas_dim-900.fits')

# 显示.fits文件的内容
fits_file.info()

# 访问.fits文件的不同数据单元（HDU）
# HDU（Header Data Unit）是.fits文件中的数据块
# 如果你的.fits文件包含多个HDUs，你可以通过索引访问它们
# 第一个HDU通常是主数据单元（Primary HDU）

# 获取主数据单元（Primary HDU）的数据
primary_hdu = fits_file[0]
data = primary_hdu.data

# 关闭.fits文件
fits_file.close()

# 保存Numpy数组为.fits文件
fits.writeto('cr2096.fits', data, overwrite=True)
