#!/usr/bin/env python3
import pycvm

def test_constants():
    """
    Verify the CVM model constants.
    """
    c = pycvm.Constants()
    assert c.get_number_of_layers() == 3, 'number of layers is wrong'
    assert c.get_utm_zone() == 10, 'utm zone is wrong' 
   
    assert c.get_number_of_grid_points_in_x(pycvm.LayerIdentifier.top) == 3271, 'nx top wrong'
    assert c.get_number_of_grid_points_in_x(pycvm.LayerIdentifier.middle) == 2181, 'nx middle wrong'
    assert c.get_number_of_grid_points_in_x(pycvm.LayerIdentifier.bottom) == 727, 'nx bottom wrong'

    assert c.get_number_of_grid_points_in_y(pycvm.LayerIdentifier.top) == 5367, 'ny top wrong'
    assert c.get_number_of_grid_points_in_y(pycvm.LayerIdentifier.middle) == 3578, 'ny middle wrong'
    assert c.get_number_of_grid_points_in_y(pycvm.LayerIdentifier.bottom) == 1193, 'ny bottom wrong'

    assert c.get_number_of_grid_points_in_z(pycvm.LayerIdentifier.top) == 13, 'nz top wrong'
    assert c.get_number_of_grid_points_in_z(pycvm.LayerIdentifier.middle) == 29, 'nz middle wrong'
    assert c.get_number_of_grid_points_in_z(pycvm.LayerIdentifier.bottom) == 55, 'nz bottom wrong'

    assert abs(c.get_grid_spacing_in_x(pycvm.LayerIdentifier.top) - 200) < 1.e-13, 'dx top wrong'
    assert abs(c.get_grid_spacing_in_x(pycvm.LayerIdentifier.middle) - 300) < 1.e-13, 'dx middle wrong'
    assert abs(c.get_grid_spacing_in_x(pycvm.LayerIdentifier.bottom) - 900) < 1.e-13, 'dx bottom wrong'

    assert abs(c.get_grid_spacing_in_y(pycvm.LayerIdentifier.top) - 200) < 1.e-13, 'dy top wrong'
    assert abs(c.get_grid_spacing_in_y(pycvm.LayerIdentifier.middle) - 300) < 1.e-13, 'dy middle wrong'
    assert abs(c.get_grid_spacing_in_y(pycvm.LayerIdentifier.bottom) - 900) < 1.e-13, 'dy bottom wrong'

    assert abs(c.get_grid_spacing_in_z(pycvm.LayerIdentifier.top) - 100) < 1.e-13, 'dz top wrong'
    assert abs(c.get_grid_spacing_in_z(pycvm.LayerIdentifier.middle) - 300) < 1.e-13, 'dz middle wrong'
    assert abs(c.get_grid_spacing_in_z(pycvm.LayerIdentifier.bottom) - 900) < 1.e-13, 'dz bottom wrong'

    """ 
    auto mLCM = static_cast<double> (std::lcm(std::lcm(200, 300), 900));
    EXPECT_NEAR(c.getLeastCommonMultipleOfGridSpacingsInXAndY(), mLCM, 1.e-12);
    """
    print("Passed constants test")

def test_geodetic():
    geo = pycvm.Geodetic()

    tol = 10. # 10 meters is more than close enough
    seattle = (47.6062, -122.3321)
    seattleUTM = (550200, 5272748.59)
    xy = geo.latitude_longitude_to_utm(seattle)
    assert abs(xy[0] - seattleUTM[0]) < tol, 'x wrong'
    assert abs(xy[1] - seattleUTM[1]) < tol, 'y wrong'
    assert abs(xy[0] - 550200.21335878095) < 1.e-7, 'x high precision is wrong' # from fortran code
    assert abs(xy[1] - 5272751.6364032440) < 1.e-7, 'y high precision is wrong'
    ll = geo.utm_to_latitude_longitude(xy)
    assert abs(ll[0] - seattle[0]) < 1.e-4, 'lat wrong'
    assert abs(ll[1] - seattle[1]) < 1.e-4, 'lon wrong'
    assert abs(ll[0] - 47.606227487491026) < 1.e-7, 'lat high precision is wrong' # From Fortran code
    assert abs(ll[1] - -122.33209965126906) < 1.e-7, 'lon high precision is wrong'
    print("Passed geodetic test")


if __name__ == "__main__":
    test_constants()
    test_geodetic()
