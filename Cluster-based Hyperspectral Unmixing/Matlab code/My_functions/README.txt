##########################
#                        #
# My_functions - Content #
# by A.Boselli           #
##########################


%%----------------------------------------------------------------------------------------

check_abundances_sum.m - Plot the map of the pixelwise sum of the abundance

It relies on plot_cube_3D.m


%%----------------------------------------------------------------------------------------

check_ANC_violations.m - Plot the map of ANC violations of the abundances cube

It relies on plot_cube_3D.m


%%----------------------------------------------------------------------------------------

check_ASC_violations.m - Plot the map of ASC violations of the abundances cube

It relies on plot_cube_3D.m


%%----------------------------------------------------------------------------------------

number_active_endmembers.m - Plot the map of the number of active endmembers per pixel

It relies on plot_cube_3D.m


%%----------------------------------------------------------------------------------------

plot_cube_3D.m - Plot all the planes of a datacube

Each plane is plotted with a number because its first use is for maps dependent on the endmember


%%----------------------------------------------------------------------------------------

plot_total_abundances.m - Barplot of the total abundances of each endmember


%%----------------------------------------------------------------------------------------

reshape_into_2D.m - Reshape a datacube into a 2D tensor (matrix)


%%----------------------------------------------------------------------------------------

reshape_into_3D.m - Reshape an image cast as a matrix into a datacube

reshape_into_2D.m and reshape_into_3D.m are meant as inverse operations on the hyperspectral image and on the datacubes in general