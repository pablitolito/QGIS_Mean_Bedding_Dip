# QGIS_Mean_Bedding_Dip
###Calculates the mean dip values using Fisher statistics and the Bingham Axial Distribution.

Each average is computed for the points contained within each polygon, provided there is more than one data point.
Each Mean_bedding point is located at the average position of the points used to calculate the mean

Input:
    layer Areas: a vector polygon layer called 'Areas', loaded in the project, with the field 'Name' indicating the name of the polygon
    layer Bedding: a vector point layer called 'Bedding', loaded in the project, with the fields 'DipDir' and 'Dip', containing dip data
Output:
    layer Mean_Bedding: a vector point layer, with the average dip values for each polygon, the pole values of the plane, and the calculated statistical parameters
    layer Unused_Points: a vector point layer, containing those points that have not been used to calculate the mean (i.e., points not inside any polygons, or points that are alone in a polygon)
    Mean_Bedding.csv file saved in the project folder, containing the same data as the Mean_bedding layer
    woodcock_diagram editable file
