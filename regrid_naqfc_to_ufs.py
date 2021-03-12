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
    #North
    import xarray as xr
    y = xr.open_dataset(fname).rename({'grid_xt':'x','grid_yt':'y','grid_lont':'lon','grid_latt':'lat'})
    lon = y.lon[0,:]
    lat = y.lat[0,:]
    return lon,lat
    
def get_top_latlon(fname='lbcs/grid_spec.nc'):
    #South
    y = xr.open_dataset(fname).rename({'grid_xt':'x','grid_yt':'y','grid_lont':'lon','grid_latt':'lat'})
    lon = y.lon[-1,:]
    lat = y.lat[-1,:]
    return lon,lat

def get_right_latlon(fname='lbcs/grid_spec.nc'):
    #West
    y = xr.open_dataset(fname).rename({'grid_xt':'x','grid_yt':'y','grid_lont':'lon','grid_latt':'lat'})
    lon = y.lon[:,-1]
    lat = y.lat[:,-1]
    return lon,lat

def get_left_latlon(fname='lbcs/grid_spec.nc'):
    #East
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
    pres = [2.0000000e-01, 6.4247000e-01, 1.3779001e+00, 2.2195797e+00, 3.1826599e+00,
       4.2843404e+00, 5.5442395e+00, 6.9845676e+00, 8.6305780e+00, 1.0510799e+01,
       1.2657519e+01, 1.5107114e+01, 1.7900507e+01, 2.1083654e+01, 2.4707880e+01,
       2.8830378e+01, 3.3514606e+01, 3.8830524e+01, 4.4854927e+01, 5.1671455e+01,
       5.9370506e+01, 6.8048744e+01, 7.7806702e+01, 8.8735672e+01, 1.0092328e+02,
       1.1445748e+02, 1.2942398e+02, 1.4590292e+02, 1.6396469e+02, 1.8366566e+02,
       2.0504283e+02, 2.2810867e+02, 2.5284596e+02, 2.7920679e+02, 3.0711218e+02,
       3.3644269e+02, 3.6703751e+02, 3.9869550e+02, 4.3117953e+02, 4.6422281e+02,
       4.9753726e+02, 5.3082397e+02, 5.6378381e+02, 5.9612842e+02, 6.2759058e+02,
       6.5793237e+02, 6.8695209e+02, 7.1448804e+02, 7.4041998e+02, 7.6466888e+02,
       7.8719440e+02, 8.0799060e+02, 8.2708185e+02, 8.4451685e+02, 8.6036389e+02,
       8.7470532e+02, 8.8763312e+02, 8.9924530e+02, 9.0964209e+02, 9.1892334e+02,
       9.2718591e+02, 9.3452368e+02, 9.4102484e+02, 9.4677258e+02, 9.5184485e+02] 
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
    pres = [2.0000000e-01, 6.4247000e-01, 1.3779001e+00, 2.2195797e+00, 3.1826599e+00,
       4.2843404e+00, 5.5442395e+00, 6.9845676e+00, 8.6305780e+00, 1.0510799e+01,
       1.2657519e+01, 1.5107114e+01, 1.7900507e+01, 2.1083654e+01, 2.4707880e+01,
       2.8830378e+01, 3.3514606e+01, 3.8830524e+01, 4.4854927e+01, 5.1671455e+01,
       5.9370506e+01, 6.8048744e+01, 7.7807938e+01, 8.8750145e+01, 1.0097821e+02,
       1.1459538e+02, 1.2970251e+02, 1.4639444e+02, 1.6475577e+02, 1.8485583e+02,
       2.0674359e+02, 2.3044174e+02, 2.5594101e+02, 2.8319666e+02, 3.1212476e+02,
       3.4259805e+02, 3.7444452e+02, 4.0744962e+02, 4.4136017e+02, 4.7589117e+02,
       5.1073523e+02, 5.4557318e+02, 5.8008557e+02, 6.1396460e+02, 6.4692474e+02,
       6.7871234e+02, 7.0911182e+02, 7.3795056e+02, 7.6510022e+02, 7.9047650e+02,
       8.1403607e+02, 8.3577295e+02, 8.5571289e+02, 8.7390784e+02, 8.9043103e+02,
       9.0536993e+02, 9.1882330e+02, 9.3089532e+02, 9.4169312e+02, 9.5132330e+02,
       9.5988953e+02, 9.6749170e+02, 9.7422412e+02, 9.8017487e+02, 9.8542584e+02] 
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
    pres = [2.00000003e-01, 6.42470419e-01, 1.37790024e+00, 2.21957755e+00,
       3.18265724e+00, 4.28433704e+00, 5.54424620e+00, 6.98456764e+00,
       8.63057613e+00, 1.05107870e+01, 1.26575060e+01, 1.51071177e+01,
       1.79005146e+01, 2.10836430e+01, 2.47078533e+01, 2.88303471e+01,
       3.35145912e+01, 3.88305359e+01, 4.48549805e+01, 5.16714973e+01,
       5.93705521e+01, 6.80487442e+01, 7.78088760e+01, 8.87609177e+01,
       1.01019051e+02, 1.14697929e+02, 1.29909576e+02, 1.46759918e+02,
       1.65343887e+02, 1.85740768e+02, 2.08008179e+02, 2.32176407e+02,
       2.58242310e+02, 2.86163208e+02, 3.15851776e+02, 3.47174652e+02,
       3.79951874e+02, 4.13958557e+02, 4.48929779e+02, 4.84566986e+02,
       5.20548401e+02, 5.56539673e+02, 5.92206421e+02, 6.27226257e+02,
       6.61300415e+02, 6.94162964e+02, 7.25588257e+02, 7.55395752e+02,
       7.83450867e+02, 8.09665283e+02, 8.33993835e+02, 8.56429871e+02,
       8.77000916e+02, 8.95761108e+02, 9.12786804e+02, 9.28170105e+02,
       9.42013977e+02, 9.54427979e+02, 9.65524048e+02, 9.75413696e+02,
       9.84205688e+02, 9.92004456e+02, 9.98908752e+02, 1.00501031e+03,
       1.01039453e+03]
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
    pres = [2.00000003e-01, 6.42470419e-01, 1.37790024e+00, 2.21957755e+00,
       3.18265724e+00, 4.28433704e+00, 5.54424620e+00, 6.98456764e+00,
       8.63057613e+00, 1.05107870e+01, 1.26575060e+01, 1.51071177e+01,
       1.79005146e+01, 2.10836430e+01, 2.47078533e+01, 2.88303471e+01,
       3.35145912e+01, 3.88305359e+01, 4.48549805e+01, 5.16714973e+01,
       5.93705521e+01, 6.80487442e+01, 7.78090057e+01, 8.87633896e+01,
       1.01028435e+02, 1.14721474e+02, 1.29957169e+02, 1.46843872e+02,
       1.65478989e+02, 1.85944046e+02, 2.08298645e+02, 2.32574860e+02,
       2.58770905e+02, 2.86844635e+02, 3.16707886e+02, 3.48225952e+02,
       3.81216919e+02, 4.15453644e+02, 4.50668488e+02, 4.86559845e+02,
       5.22802490e+02, 5.59058716e+02, 5.94990601e+02, 6.30272522e+02,
       6.64602478e+02, 6.97711975e+02, 7.29372864e+02, 7.59402893e+02,
       7.87665955e+02, 8.14072998e+02, 8.38578064e+02, 8.61174866e+02,
       8.81890686e+02, 9.00780640e+02, 9.17921875e+02, 9.33407166e+02,
       9.47340942e+02, 9.59833435e+02, 9.70997986e+02, 9.80947266e+02,
       9.89791016e+02, 9.97635071e+02, 1.00457874e+03, 1.01071515e+03,
       1.01612976e+03] 
    lbcs['pres'] = (('z'),pres)
#    print(lbcs)
#    print(lon)
#    print(lat)

    lbcs = lbcs.expand_dims('x')
    lbcs['lon'] = lon.expand_dims('x')
    lbcs['lat'] = lat.expand_dims('x')
    return lbcs.set_coords(['lon','lat','pres']).transpose('z','y','x')


def get_fv3_rough_pres():
#    This is the top pressure level
    pres = [2.0000000e-01, 6.4247000e-01, 1.3779001e+00, 2.2195797e+00, 3.1826599e+00,
       4.2843404e+00, 5.5442395e+00, 6.9845676e+00, 8.6305780e+00, 1.0510799e+01,
       1.2657519e+01, 1.5107114e+01, 1.7900507e+01, 2.1083654e+01, 2.4707880e+01,
       2.8830378e+01, 3.3514606e+01, 3.8830524e+01, 4.4854927e+01, 5.1671455e+01,
       5.9370506e+01, 6.8048744e+01, 7.7807938e+01, 8.8750145e+01, 1.0097821e+02,
       1.1459538e+02, 1.2970251e+02, 1.4639444e+02, 1.6475577e+02, 1.8485583e+02,
       2.0674359e+02, 2.3044174e+02, 2.5594101e+02, 2.8319666e+02, 3.1212476e+02,
       3.4259805e+02, 3.7444452e+02, 4.0744962e+02, 4.4136017e+02, 4.7589117e+02,
       5.1073523e+02, 5.4557318e+02, 5.8008557e+02, 6.1396460e+02, 6.4692474e+02,
       6.7871234e+02, 7.0911182e+02, 7.3795056e+02, 7.6510022e+02, 7.9047650e+02,
       8.1403607e+02, 8.3577295e+02, 8.5571289e+02, 8.7390784e+02, 8.9043103e+02,
       9.0536993e+02, 9.1882330e+02, 9.3089532e+02, 9.4169312e+02, 9.5132330e+02,
       9.5988953e+02, 9.6749170e+02, 9.7422412e+02, 9.8017487e+02, 9.8542584e+02]
    return pres # xr.DataArray((('lev'),pres))
