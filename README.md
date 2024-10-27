# QGIS Mean Bedding (Dip Direction & Dip)

### Description
Calculates the mean dip values using Fisher statistics and the Bingham Axial Distribution.

Each average is computed for the points contained within each polygon, provided there is more than one data point. Each Mean_bedding point is located at the average position of the points used to calculate the mean.

### Input
- **Layer Areas**: A vector polygon layer called **'Areas'**, loaded in the project, with the field **'Name'** indicating the name of the polygon.
- **Layer Bedding**: A vector point layer called **'Bedding'**, loaded in the project, with the fields **'DipDir'** and **'Dip'**, containing dip data.

### Output
- **Layer Mean_Bedding**: A vector point layer with the average dip values for each polygon, the pole values of the plane, and the calculated statistical parameters.
- **Layer Unused_Points**: A vector point layer containing those points that have not been used to calculate the mean (i.e., points not inside any polygons, or points that are alone in a polygon).
- **Mean_Bedding.csv**: A file saved in the project folder, containing the same data as the Mean_bedding layer.
- **Woodcock diagram**: An editable file.
