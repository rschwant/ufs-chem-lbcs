# coding: utf-8
import xarray as xr
import monet as m


def fix_all(fname,cmaq_fname):
    o = open_cmaq_bl(lbcf=cmaq_fname)
    print('interpolating bottom')
    bottom = fix_lbcs_bottom(o,fname)
    bottom = bottom.drop([n for n in bottom.data_vars if 'bottom' not in n])
    print('interpolating top')
    top = fix_lbcs_top(o,fname)
    top = top.drop([n for n in top.data_vars if 'top' not in n])
    print('interpolating right')
    right = fix_lbcs_right(o,fname)
    right = right.drop([n for n in right.data_vars if 'right' not in n])
    print('interpolating left')
    left = fix_lbcs_left(o,fname)
    left = left.drop([n for n in left.data_vars if 'left' not in n])
    
    a = xr.merge([bottom,top,left,right])
    
    return a

    
#    return right,left
def fix_lbcs_top(cmaq_lbcs,fname):
    orig = xr.open_dataset(fname)
    lbcs_top = open_fv3_lbcs_for_top(fname=fname)
    lbcs = interp_top(lbcs_top,cmaq_lbcs)
    top_vars = [n for n in lbcs.data_vars if 'top' in n]
    for h in range(4):
    #h = 0
    #n = 'o3_top'
        for n in top_vars:
            newvar = orig[n].data
            newvar[:,h,4:-4] = lbcs[n].squeeze().data
            newvar[:,h,0] = lbcs[n].squeeze().data[:,0]
            newvar[:,h,1] = lbcs[n].squeeze().data[:,0]
            newvar[:,h,2] = lbcs[n].squeeze().data[:,0]
            newvar[:,h,3] = lbcs[n].squeeze().data[:,0]
            newvar[:,h,-1] = lbcs[n].squeeze().data[:,-1]
            newvar[:,h,-2] = lbcs[n].squeeze().data[:,-1]
            newvar[:,h,-3] = lbcs[n].squeeze().data[:,-1]
            newvar[:,h,-4] = lbcs[n].squeeze().data[:,-1]
            orig[n].data = newvar
            orig[n].data = orig[n].where(orig[n] > 0,0)
    dropvars = [n for n    in lbcs.data_vars if 'top' not in n]

    return orig.drop(dropvars)

def fix_lbcs_left(cmaq_lbcs,fname):
    orig = xr.open_dataset(fname)
    lbcs_left = open_fv3_lbcs_for_left(fname=fname)
    lbcs = interp_left(lbcs_left,cmaq_lbcs)
    left_vars = [n for n in lbcs.data_vars if 'left' in n]
    for h in range(4):
    #h = 0
    #n = 'o3_left'
        for n in left_vars:
            newvar = orig[n].data
            newvar[:,:,h] = lbcs[n].squeeze().data
            orig[n].data = newvar
            orig[n].data = orig[n].where(orig[n] > 0,0)
    dropvars = [n for n    in lbcs.data_vars if 'left' not in n]

    return orig.drop(dropvars)

def fix_lbcs_right(cmaq_lbcs,fname):
    orig = xr.open_dataset(fname)
    lbcs_right = open_fv3_lbcs_for_right(fname=fname)
    lbcs = interp_right(lbcs_right,cmaq_lbcs)
    right_vars = [n for n in lbcs.data_vars if 'right' in n]
    for h in range(4):
    #h = 0
    #n = 'o3_right'
        for n in right_vars:
            newvar = orig[n].data
            newvar[:,:,h] = lbcs[n].squeeze().data
            orig[n].data = newvar
            orig[n].data = orig[n].where(orig[n] > 0,0)
    dropvars = [n for n    in lbcs.data_vars if 'right' not in n]

    return orig.drop(dropvars)

def fix_lbcs_bottom(cmaq_lbcs,fname):
    orig = xr.open_dataset(fname)
    lbcs_bottom = open_fv3_lbcs_for_bottom(fname=fname)
    lbcs = interp_bottom(lbcs_bottom,cmaq_lbcs)
    bottom_vars = [n for n in lbcs.data_vars if 'bottom' in n]
    for h in range(4):
    #h = 0
    #n = 'o3_bottom'
        for n in bottom_vars:
            newvar = orig[n].data
            newvar[:,h,4:-4] = lbcs[n].squeeze().data
            newvar[:,h,0] = lbcs[n].squeeze().data[:,0]
            newvar[:,h,1] = lbcs[n].squeeze().data[:,0]
            newvar[:,h,2] = lbcs[n].squeeze().data[:,0]
            newvar[:,h,3] = lbcs[n].squeeze().data[:,0]
            newvar[:,h,-1] = lbcs[n].squeeze().data[:,-1]
            newvar[:,h,-2] = lbcs[n].squeeze().data[:,-1]
            newvar[:,h,-3] = lbcs[n].squeeze().data[:,-1]
            newvar[:,h,-4] = lbcs[n].squeeze().data[:,-1]
            orig[n].data = newvar
            orig[n].data = orig[n].where(orig[n] > 0,0)
    dropvars = [n for n in lbcs.data_vars if 'bottom' not in n]
    
    return orig.drop(dropvars)

    
    
def interp_cmaq_to_fv3_pres(cmaq_pres,fv3_pres,cmaq_val):
    from scipy.interpolate import interp1d
    lower,upper = cmaq_val[0],cmaq_val[-1]
    f1 = interp1d(cmaq_pres,cmaq_val,fill_value=(upper,lower),bounds_error=False,kind='cubic')
    ynew = f1(fv3_pres)
    return ynew
    
    
def interp_bottom(lbcs_bottom,cmaq_lbcs):
    hmapped = lbcs_bottom.o3_bottom.monet.remap_nearest(cmaq_lbcs,radius_of_influence=1000000).compute()
    for n in hmapped.data_vars:
        cmaq_var = hmapped[n]
        fv3_var = lbcs_bottom[n + '_bottom'].copy()
        for i in range(len(hmapped.x)):
            ynew = interp_cmaq_to_fv3_pres(cmaq_var.pres.values,fv3_var.pres.squeeze().values,cmaq_var.isel(y=0,x=i).squeeze())
    #        print(ynew)
            fv3_var[:,0,i] = ynew
        lbcs_bottom[n + '_bottom'] = fv3_var
    return lbcs_bottom

def interp_top(lbcs_bottom,cmaq_lbcs):
    hmapped = lbcs_bottom.o3_top.monet.remap_nearest(cmaq_lbcs,radius_of_influence=1000000).compute()
    for n in hmapped.data_vars:
        cmaq_var = hmapped[n]
        fv3_var = lbcs_bottom[n + '_top'].copy()
        for i in range(len(hmapped.x)):
            ynew = interp_cmaq_to_fv3_pres(cmaq_var.pres.values,fv3_var.pres.squeeze().values,cmaq_var.isel(y=0,x=i).squeeze())
    #        print(ynew)
            fv3_var[:,0,i] = ynew
        lbcs_bottom[n + '_top'] = fv3_var
    return lbcs_bottom
 
def interp_left(lbcs_bottom,cmaq_lbcs):
    hmapped = lbcs_bottom.o3_left.monet.remap_nearest(cmaq_lbcs,radius_of_influence=1000000).compute()
    for n in hmapped.data_vars:
        cmaq_var = hmapped[n]
        fv3_var = lbcs_bottom[n + '_left'].copy()
        for i in range(len(hmapped.y)):
            ynew = interp_cmaq_to_fv3_pres(cmaq_var.pres.values,fv3_var.pres.squeeze().values,cmaq_var.isel(y=i,x=0).squeeze())
    #        print(ynew)
            fv3_var[:,i,0] = ynew
        lbcs_bottom[n + '_left'] = fv3_var
    return lbcs_bottom

def interp_right(lbcs_bottom,cmaq_lbcs):
    hmapped = lbcs_bottom.o3_right.monet.remap_nearest(cmaq_lbcs,radius_of_influence=1000000).compute()
    for n in hmapped.data_vars:
        cmaq_var = hmapped[n]
        fv3_var = lbcs_bottom[n + '_right'].copy()
        for i in range(len(hmapped.y)):
            ynew = interp_cmaq_to_fv3_pres(cmaq_var.pres.values,fv3_var.pres.squeeze().values,cmaq_var.isel(y=i,x=0).squeeze())
    #        print(ynew)
            fv3_var[:,i,0] = ynew
        lbcs_bottom[n + '_right'] = fv3_var
    return lbcs_bottom   
def get_bottom_latlon(fname='lbcs/grid_spec.nc'):
    import xarray as xr
    y = xr.open_dataset(fname).rename({'grid_xt':'x','grid_yt':'y','grid_lont':'lon','grid_latt':'lat'})
    lon = y.lon[-1,:]
    lat = y.lat[-1,:]
    return lon,lat
    
def get_top_latlon(fname='lbcs/grid_spec.nc'):
    y = xr.open_dataset(fname).rename({'grid_xt':'x','grid_yt':'y','grid_lont':'lon','grid_latt':'lat'})
    lon = y.lon[0,:]
    lat = y.lat[0,:]
    return lon,lat

def get_right_latlon(fname='lbcs/grid_spec.nc'):
    y = xr.open_dataset(fname).rename({'grid_xt':'x','grid_yt':'y','grid_lont':'lon','grid_latt':'lat'})
    lon = y.lon[:,-1]
    lat = y.lat[:,-1]
    return lon,lat

def get_left_latlon(fname='lbcs/grid_spec.nc'):
    y = xr.open_dataset(fname).rename({'grid_xt':'x','grid_yt':'y','grid_lont':'lon','grid_latt':'lat'})
    lon = y.lon[:,0]
    lat = y.lat[:,0]
    return lon,lat

def open_cmaq_bl(lbcf='aqm_conus_12km_geos_200608_static_35L.ncf',lbcgrid='aqm.t00z.grdbdy2d.ncf',lbcmet='aqm.t00z.metbdy3d.ncf'):
    o = xr.open_dataset(lbcf,drop_variables=['TFLAG']).squeeze()
    om = xr.open_dataset(lbcmet).isel(TSTEP=0).mean(dim='PERIM').squeeze()
    og = xr.open_dataset(lbcgrid).squeeze()
    o['PRES'] = om['PRES']
    o['LON'] = og['LON']
    o['LAT'] = og['LAT']
    rename_dict = {}
    for n in o.data_vars:
        rename_dict[n] = n.lower()
    o = o.rename(rename_dict)
    o = o.rename({'PERIM':'x','LAY':'z'})
    o = o.expand_dims('y').transpose('z','y','x')
    o['pres'] = o.pres.squeeze()/100
    return o.set_coords(['lon','lat','pres'])

def open_fv3_lbcs_for_bottom(fname='lbcs/gfs_bndy_chem_08.tile7.000.nc',grid_fname='lbcs/grid_spec.nc'):
    import xarray as xr
    lon,lat=get_bottom_latlon(grid_fname)
    lbcs = xr.open_dataset(fname).isel(halo=0)
    lbcs = lbcs.rename({'lon':'x','lat':'y','lev':'z'})
    drop_vars = [n for n in lbcs.data_vars if 'bottom' not in n]
    lbcs = lbcs.drop(drop_vars)
#    lbcs['lon'] = lon
#    lbcs['lat'] = lat
 #   return lbcs.set_coords(['lon','lat'])
    lbcs = lbcs.isel(x=slice(4,len(lbcs.x) - 4))
    pres = [1.993046e-01, 5.608634e+00, 1.173876e+01, 1.867534e+01, 2.651118e+01,
       3.534573e+01, 4.528432e+01, 5.643637e+01, 6.891479e+01, 8.283226e+01,
       9.829910e+01, 1.154193e+02, 1.342857e+02, 1.549752e+02, 1.775427e+02,
       2.020157e+02, 2.283876e+02, 2.566128e+02, 2.866021e+02, 3.182203e+02,
       3.512849e+02, 3.855683e+02, 4.208025e+02, 4.566861e+02, 4.928940e+02,
       5.290892e+02, 5.649343e+02, 6.001047e+02, 6.342996e+02, 6.672516e+02,
       6.987342e+02, 7.285657e+02, 7.566116e+02, 7.827839e+02, 8.070405e+02,
       8.293792e+02, 8.498318e+02, 8.684576e+02, 8.853377e+02, 9.005690e+02,
       9.142595e+02, 9.265233e+02, 9.374772e+02, 9.472375e+02, 9.559158e+02,
       9.636148e+02, 9.704313e+02, 9.764557e+02, 9.817717e+02, 9.864559e+02,
       9.905784e+02, 9.942026e+02, 9.973857e+02, 1.000179e+03, 1.002629e+03,
       1.004775e+03, 1.006656e+03, 1.008302e+03, 1.009742e+03, 1.011002e+03,
       1.012104e+03, 1.013067e+03, 1.013908e+03, 1.014644e+03,101508.617188/100]
    lbcs['pres'] = (('z'),pres)
    lbcs = lbcs.expand_dims('y')
    lbcs['lon'] = lon.expand_dims('y')
    lbcs['lat'] = lat.expand_dims('y')
    return lbcs.set_coords(['lon','lat','pres']).transpose('z','y','x')
    
def open_fv3_lbcs_for_top(fname='lbcs/gfs_bndy_chem_08.tile7.000.nc',grid_fname='lbcs/grid_spec.nc'):
    import xarray as xr
    lon,lat=get_top_latlon(grid_fname)
    lbcs = xr.open_dataset(fname).isel(halo=0)
    lbcs = lbcs.rename({'lon':'x','lat':'y','lev':'z'})
    drop_vars = [n for n in lbcs.data_vars if 'top' not in n]
    lbcs = lbcs.drop(drop_vars)
#    lbcs['lon'] = lon
#    lbcs['lat'] = lat
 #   return lbcs.set_coords(['lon','lat'])
    lbcs = lbcs.isel(x=slice(4,len(lbcs.x) - 4))
    pres = [1.993046e-01, 5.608634e+00, 1.173876e+01, 1.867534e+01, 2.651118e+01,
       3.534573e+01, 4.528432e+01, 5.643637e+01, 6.891479e+01, 8.283226e+01,
       9.829910e+01, 1.154193e+02, 1.342857e+02, 1.549752e+02, 1.775427e+02,
       2.020157e+02, 2.283876e+02, 2.566128e+02, 2.866021e+02, 3.182203e+02,
       3.512849e+02, 3.855683e+02, 4.208025e+02, 4.566861e+02, 4.928940e+02,
       5.290892e+02, 5.649343e+02, 6.001047e+02, 6.342996e+02, 6.672516e+02,
       6.987342e+02, 7.285657e+02, 7.566116e+02, 7.827839e+02, 8.070405e+02,
       8.293792e+02, 8.498318e+02, 8.684576e+02, 8.853377e+02, 9.005690e+02,
       9.142595e+02, 9.265233e+02, 9.374772e+02, 9.472375e+02, 9.559158e+02,
       9.636148e+02, 9.704313e+02, 9.764557e+02, 9.817717e+02, 9.864559e+02,
       9.905784e+02, 9.942026e+02, 9.973857e+02, 1.000179e+03, 1.002629e+03,
       1.004775e+03, 1.006656e+03, 1.008302e+03, 1.009742e+03, 1.011002e+03,
       1.012104e+03, 1.013067e+03, 1.013908e+03, 1.014644e+03,101508.617188/100]
    lbcs['pres'] = (('z'),pres)
    lbcs = lbcs.expand_dims('y')
    lbcs['lon'] = lon.expand_dims('y')
    lbcs['lat'] = lat.expand_dims('y')
    return lbcs.set_coords(['lon','lat','pres']).transpose('z','y','x')

def open_fv3_lbcs_for_left(fname='lbcs/gfs_bndy_chem_08.tile7.000.nc',grid_fname='lbcs/grid_spec.nc'):
    import xarray as xr
    lon,lat=get_left_latlon(grid_fname)
    lbcs = xr.open_dataset(fname).isel(halo=0)
    lbcs = lbcs.rename({'lon':'x','lat':'y','lev':'z'})
    drop_vars = [n for n in lbcs.data_vars if 'left' not in n]
    lbcs = lbcs.drop(drop_vars)
#    lbcs['lon'] = lon
#    lbcs['lat'] = lat
 #   return lbcs.set_coords(['lon','lat'])
#    lbcs = lbcs.isel(x=slice(4,len(lbcs.x) - 4))
    pres = [1.993046e-01, 5.608634e+00, 1.173876e+01, 1.867534e+01, 2.651118e+01,
       3.534573e+01, 4.528432e+01, 5.643637e+01, 6.891479e+01, 8.283226e+01,
       9.829910e+01, 1.154193e+02, 1.342857e+02, 1.549752e+02, 1.775427e+02,
       2.020157e+02, 2.283876e+02, 2.566128e+02, 2.866021e+02, 3.182203e+02,
       3.512849e+02, 3.855683e+02, 4.208025e+02, 4.566861e+02, 4.928940e+02,
       5.290892e+02, 5.649343e+02, 6.001047e+02, 6.342996e+02, 6.672516e+02,
       6.987342e+02, 7.285657e+02, 7.566116e+02, 7.827839e+02, 8.070405e+02,
       8.293792e+02, 8.498318e+02, 8.684576e+02, 8.853377e+02, 9.005690e+02,
       9.142595e+02, 9.265233e+02, 9.374772e+02, 9.472375e+02, 9.559158e+02,
       9.636148e+02, 9.704313e+02, 9.764557e+02, 9.817717e+02, 9.864559e+02,
       9.905784e+02, 9.942026e+02, 9.973857e+02, 1.000179e+03, 1.002629e+03,
       1.004775e+03, 1.006656e+03, 1.008302e+03, 1.009742e+03, 1.011002e+03,
       1.012104e+03, 1.013067e+03, 1.013908e+03, 1.014644e+03,101508.617188/100]
    lbcs['pres'] = (('z'),pres)
#    print(lbcs)
    lbcs = lbcs.expand_dims('x')
    lbcs['lon'] = lon.expand_dims('x')
    lbcs['lat'] = lat.expand_dims('x')
    return lbcs.set_coords(['lon','lat','pres']).transpose('z','y','x')

def open_fv3_lbcs_for_right(fname='lbcs/gfs_bndy_chem_08.tile7.000.nc',grid_fname='lbcs/grid_spec.nc'):
    import xarray as xr
    lon,lat=get_right_latlon(grid_fname)
    lbcs = xr.open_dataset(fname).isel(halo=0)
    lbcs = lbcs.rename({'lon':'x','lat':'y','lev':'z'})
    drop_vars = [n for n in lbcs.data_vars if 'right' not in n]
    lbcs = lbcs.drop(drop_vars)
#    lbcs['lon'] = lon
#    lbcs['lat'] = lat
 #   return lbcs.set_coords(['lon','lat'])
#    lbcs = lbcs.isel(x=slice(4,len(lbcs.x) - 4))
    pres = [1.993046e-01, 5.608634e+00, 1.173876e+01, 1.867534e+01, 2.651118e+01,
       3.534573e+01, 4.528432e+01, 5.643637e+01, 6.891479e+01, 8.283226e+01,
       9.829910e+01, 1.154193e+02, 1.342857e+02, 1.549752e+02, 1.775427e+02,
       2.020157e+02, 2.283876e+02, 2.566128e+02, 2.866021e+02, 3.182203e+02,
       3.512849e+02, 3.855683e+02, 4.208025e+02, 4.566861e+02, 4.928940e+02,
       5.290892e+02, 5.649343e+02, 6.001047e+02, 6.342996e+02, 6.672516e+02,
       6.987342e+02, 7.285657e+02, 7.566116e+02, 7.827839e+02, 8.070405e+02,
       8.293792e+02, 8.498318e+02, 8.684576e+02, 8.853377e+02, 9.005690e+02,
       9.142595e+02, 9.265233e+02, 9.374772e+02, 9.472375e+02, 9.559158e+02,
       9.636148e+02, 9.704313e+02, 9.764557e+02, 9.817717e+02, 9.864559e+02,
       9.905784e+02, 9.942026e+02, 9.973857e+02, 1.000179e+03, 1.002629e+03,
       1.004775e+03, 1.006656e+03, 1.008302e+03, 1.009742e+03, 1.011002e+03,
       1.012104e+03, 1.013067e+03, 1.013908e+03, 1.014644e+03,101508.617188/100]
    lbcs['pres'] = (('z'),pres)
#    print(lbcs)
#    print(lon)
#    print(lat)

    lbcs = lbcs.expand_dims('x')
    lbcs['lon'] = lon.expand_dims('x')
    lbcs['lat'] = lat.expand_dims('x')
    return lbcs.set_coords(['lon','lat','pres']).transpose('z','y','x')


def get_fv3_rough_pres():
    pres = [1.993046e-01, 5.608634e+00, 1.173876e+01, 1.867534e+01, 2.651118e+01,
            3.534573e+01, 4.528432e+01, 5.643637e+01, 6.891479e+01, 8.283226e+01,
            9.829910e+01, 1.154193e+02, 1.342857e+02, 1.549752e+02, 1.775427e+02,
            2.020157e+02, 2.283876e+02, 2.566128e+02, 2.866021e+02, 3.182203e+02,
            3.512849e+02, 3.855683e+02, 4.208025e+02, 4.566861e+02, 4.928940e+02,
            5.290892e+02, 5.649343e+02, 6.001047e+02, 6.342996e+02, 6.672516e+02,
            6.987342e+02, 7.285657e+02, 7.566116e+02, 7.827839e+02, 8.070405e+02,
            8.293792e+02, 8.498318e+02, 8.684576e+02, 8.853377e+02, 9.005690e+02,
            9.142595e+02, 9.265233e+02, 9.374772e+02, 9.472375e+02, 9.559158e+02,
            9.636148e+02, 9.704313e+02, 9.764557e+02, 9.817717e+02, 9.864559e+02,
            9.905784e+02, 9.942026e+02, 9.973857e+02, 1.000179e+03, 1.002629e+03,
            1.004775e+03, 1.006656e+03, 1.008302e+03, 1.009742e+03, 1.011002e+03,
            1.012104e+03, 1.013067e+03, 1.013908e+03, 1.014644e+03,101508.617188/100]
    return pres # xr.DataArray((('lev'),pres))
