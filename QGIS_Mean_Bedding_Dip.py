from qgis.core import *
from qgis.PyQt.QtCore import *
import numpy as np
from numpy.linalg import eig
import os
from qgis.core import QgsVectorFileWriter
import matplotlib.pyplot as plt
import matplotlib.cm as cm

"""
Calculates the mean dip values using Fisher statistics and the Bingham Axial Distribution.
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
"""

def fisher_mean(dir):
    '''Calculates Fisher mean'''
    fisher_mean = []
    deg = 180 / np.pi
    car = dir2car(dir)
    N = len(dir)
    X, Y, Z = 0, 0, 0
    
    for site in car:
        X += site[0]
        Y += site[1]
        Z += site[2]

    R = np.sqrt(X**2 + Y**2 + Z**2)
    k = (N - 1) / (N - R)
    
    cosAlfa = 1 - ((N - R) / R) * (20.0 ** (1.0 / (N - 1)) - 1)
    if cosAlfa < -1: cosAlfa = -1
    alpha95 = deg * np.arccos(cosAlfa)

    AzMean = deg * np.arctan2(Y, X)
    '''
    if X < 0.0:
        AzMean += 180.0
    elif Y < 0.0:
        AzMean += 360.0
    '''
    IncMean = deg * np.arcsin(Z / R)
    
    fisher_mean = [AzMean, IncMean, k, alpha95]
    return AzMean, IncMean, k, alpha95

def dir2car(dir):
    '''Converts geographic to Cartesian coordinates'''
    cart = []
    rad = np.pi / 180

    for site in dir:
        x = np.cos(rad * site[0]) * np.cos(rad * site[1]) 
        y = np.sin(rad * site[0]) * np.cos(rad * site[1]) 
        z = np.sin(rad * site[1])                 
        cartSite = [x, y, z]
        cart.append(cartSite)
    return cart


def az_inc_to_cartesian(az, inc):
    '''Function to convert azimuth and inclination to direction cosines'''
    az_rad = np.radians(az)
    inc_rad = np.radians(inc)
    x = np.cos(inc_rad) * np.sin(az_rad)
    y = np.cos(inc_rad) * np.cos(az_rad)
    z = np.sin(inc_rad)
    return np.array([x, y, z])

def bingham_orientation_matrix(az_inc_list):
    '''Calculation of the Bingham orientation matrix'''
    T = np.zeros((3, 3))  # Matriz de orientación
    for az, inc in az_inc_list:
        p = az_inc_to_cartesian(az, inc)
        T += np.outer(p, p)  # Suma de productos externos
    T /= len(az_inc_list)  # Normalizar por el número de datos
    return T

def eigen_decomposition(T):
    '''Calculation of eigenvalues and eigenvectors'''
    eigenvalues, eigenvectors = eig(T)
    order = np.argsort(eigenvalues)[::-1]  # Ordenar de mayor a menor
    eigenvalues = eigenvalues[order]
    eigenvectors = eigenvectors[:, order]
    return eigenvalues, eigenvectors

def eigenvector_to_az_inc(eigenvector):
    '''Calculate average directions of the eigenvectors in geographic coordinates'''
    x, y, z = eigenvector
    az = np.degrees(np.arctan2(x, y)) % 360
    inc = np.degrees(np.arcsin(z))
    return az, inc

def process_orientation_data(az_inc_list):
    '''Main function to process the data'''
    # Fisher mean
    az_fisher, inc_fisher, kappa, alpha95 = fisher_mean(az_inc_list)
    
    # Bingham orientation matrix
    T = bingham_orientation_matrix(az_inc_list)
    
    # Eigenvalues and eigenvectors
    eigenvalues, eigenvectors = eigen_decomposition(T)
    
    # Convert eigenvectors to geographic coordinates
    az_eigen1, inc_eigen1 = eigenvector_to_az_inc(eigenvectors[:, 0])
    az_eigen2, inc_eigen2 = eigenvector_to_az_inc(eigenvectors[:, 1])
    az_eigen3, inc_eigen3 = eigenvector_to_az_inc(eigenvectors[:, 2])
    
    return az_fisher, inc_fisher, kappa, alpha95, eigenvalues, az_eigen1, inc_eigen1, az_eigen2, inc_eigen2, az_eigen3, inc_eigen3
    
def woodcock_plot(names, eigenvalues, n):
    '''Function to draw the Woodcock plot'''
    
    eigenvalues = np.array(eigenvalues)
    n = np.array(n)

    # Calculate xy values for the plot
    x = np.abs(eigenvalues[:, 2] / eigenvalues[:, 1])  # (|E3|/|E2|)
    y = np.abs(eigenvalues[:, 1] / eigenvalues[:, 0])  # (|E2|/|E1|)

    # Create figure and adjust axes
    fig, axs = plt.subplots(1, 2, figsize=(12, 6), gridspec_kw={'width_ratios': [1, 1]})
    fig.suptitle('Woodcock Diagram', fontsize=16, fontweight='bold')
        
    # Plot with normal axes
    scatter1 = axs[0].scatter(x, y, c=n, cmap='viridis', marker='s')
    axs[0].set_xlabel('|E3/E2|')
    axs[0].set_ylabel('|E2/E1|')
    axs[0].set_title(f'n = {len(x)}')
    axs[0].axis('equal')

    # Add labels with points names
    #for i, name in enumerate(names):
    #    axs[0].annotate(name, (x[i], y[i]), textcoords="offset points", xytext=(5, 5), ha='center')

    # Draw the line x=y
    max_limit = max(np.max(x), np.max(y))
    axs[0].plot([0, max_limit], [0, max_limit], color='black', linestyle='-', linewidth=1)

    # Add labels for Girdles and Clusters
    axs[0].text(max_limit * 0.75, max_limit * 0.92, 'GIRDLES', fontsize=12, ha='center', va='center', fontweight='bold')
    axs[0].text(max_limit * 0.92, max_limit * 0.75, 'CLUSTERS', fontsize=12, ha='center', va='center', fontweight='bold')

    # Filter valid indices for the logarithmic plot
    valid_indices = np.where((eigenvalues[:, 0] > 0.0000001) & (eigenvalues[:, 1] > 0.0000001) & (eigenvalues[:, 2] > 0.0000001))[0]
    x_log = np.abs(eigenvalues[valid_indices, 2] / eigenvalues[valid_indices, 1])
    y_log = np.abs(eigenvalues[valid_indices, 1] / eigenvalues[valid_indices, 0])

    # Plot with logarithmic scales
    scatter2 = axs[1].scatter(x_log, y_log, c=n[valid_indices], cmap='viridis', marker='s')
    axs[1].set_xscale('log')
    axs[1].set_yscale('log')
    axs[1].set_xlabel('|E3/E2| (log scale)')
    axs[1].set_ylabel('|E2/E1| (log scale)')
    axs[1].set_title(f'n = {len(x_log)}, data with x&y >0')
    #axs[1].set_aspect('equal', 'box')
    axs[1].axis('equal')
    
    axs[1].plot([0, max_limit], [0, max_limit], color='black', linestyle='-', linewidth=1)

    # Add labels with points names
    #for i, name in enumerate(np.array(names)[valid_indices]):
    #    axs[1].annotate(name, (x_log[i], y_log[i]), textcoords="offset points", xytext=(5, 5), ha='center')

    # Create a color bar on the left
    cbar_ax = fig.add_axes([0.03, 0.10, 0.02, 0.7])  # Posición en la figura
    cbar = fig.colorbar(scatter1, cax=cbar_ax)
    cbar.set_label('Number of data used to calculate the mean', rotation=90, labelpad=-48)

    plt.tight_layout(rect=[0.08, 0, 1, 0.95])  # Espacio para la barra de color y el título
    
    # Save the figure
    project_path = QgsProject.instance().homePath()
    pdf_path = os.path.join(project_path, 'woodcock_diagram.pdf')
    svg_path = os.path.join(project_path, 'woodcock_diagram.svg')
    plt.savefig(pdf_path, bbox_inches='tight', pad_inches=0.1)
    plt.savefig(svg_path, bbox_inches='tight', pad_inches=0.1)
    plt.show()




def mean_bedding_main():
    # Load input layers
    areas_layer = QgsProject.instance().mapLayersByName('Areas')[0]
    bedding_layer = QgsProject.instance().mapLayersByName('Bedding')[0]

    # Check coordinate reference system
    if areas_layer.crs() != bedding_layer.crs():
        print("Atención: Las capas no están en el mismo sistema de coordenadas.")
        return

    # Create output layer for Mean_bedding
    mean_bedding_layer = QgsVectorLayer('Point?crs=EPSG:32630', 'Mean_bedding', 'memory')
    mean_bedding_provider = mean_bedding_layer.dataProvider()
    mean_bedding_provider.addAttributes([
        QgsField('Name_poly', QVariant.String),
        QgsField("X_Mean", QVariant.Double),
        QgsField("Y_Mean", QVariant.Double),
        QgsField('n', QVariant.Int),
        QgsField('DipDir_Mean', QVariant.Double),
        QgsField('Dip_Mean', QVariant.Double),
        QgsField('Pole_Az_Mean', QVariant.Double),
        QgsField('Pole_Inc_Mean', QVariant.Double),
        QgsField('kappa', QVariant.Double),
        QgsField('Alpha95', QVariant.Double),
        QgsField("E_1", QVariant.Double),
        QgsField("E_2", QVariant.Double),
        QgsField('E_3', QVariant.Double),
        QgsField('E_1_az', QVariant.Double),
        QgsField('E_1_inc', QVariant.Double),
        QgsField('E_2_az', QVariant.Double),
        QgsField('E_2_inc', QVariant.Double),
        QgsField('E_3_az', QVariant.Double),
        QgsField('E_3_inc', QVariant.Double),
    ])
    mean_bedding_layer.updateFields()

    # Create layer for unused points
    unused_points_layer = QgsVectorLayer('Point?crs=EPSG:32630', 'Unused_points', 'memory')
    unused_provider = unused_points_layer.dataProvider()
    unused_provider.addAttributes(bedding_layer.fields())
    unused_points_layer.updateFields()

    used_points = set()
    all_points_inside = []
        
    wc_names=[]
    wc_eigenvalues=[]
    wc_n=[]
    
    # Start editing to add calculated mean points
    mean_bedding_layer.startEditing()

    # Process each polygon in the 'Areas' layer
    for feature in areas_layer.getFeatures():
        geom = feature.geometry()
        points_inside = []

        # Get points inside the polygon
        for point_feature in bedding_layer.getFeatures():
            if geom.contains(point_feature.geometry()):
                dipdir = point_feature['DipDir']
                dip = point_feature['dip']
                points_inside.append([dipdir, dip, point_feature.geometry()])
                all_points_inside.append(point_feature)
                used_points.add(point_feature.id())

        # Calculate values only if there is more than one point in the polygon
        if len(points_inside) > 1:
            polygon_name = feature['Name']
            az_inc_pairs = []
            coords_sum = [0, 0]
            
            # Calculate the Azimuth and Inclination for each point
            for dipdir, dip, point_geom in points_inside:
                az = (dipdir + 180) % 360
                inc = 90 - dip
                az_inc_pairs.append([az, inc])
                coords_sum[0] += point_geom.asPoint().x()
                coords_sum[1] += point_geom.asPoint().y()

            # Get orientation parameters
            az_mean, inc_mean, kappa, alpha95, eigenvalues, az_eigen1, inc_eigen1, az_eigen2, inc_eigen2, az_eigen3, inc_eigen3 = process_orientation_data(az_inc_pairs)

            # Calculate DipDir_Mean and Dip_Mean
            dipdir_mean = (az_mean - 180) % 360
            dip_mean = 90 - inc_mean
            n = len(points_inside)

            # Calculate the centroid of the points
            centroid_x = coords_sum[0] / n
            centroid_y = coords_sum[1] / n

            # Create the mean point and assign attributes
            mean_point = QgsFeature()
            mean_point.setGeometry(QgsGeometry.fromPointXY(QgsPointXY(centroid_x, centroid_y)))
            mean_point.setAttributes([
                polygon_name,
                round(float(centroid_x), 1),
                round(float(centroid_y), 1),
                n,
                round(float(dipdir_mean), 1),
                round(float(dip_mean), 1),
                round(float(az_mean), 1),
                round(float(inc_mean), 1),
                round(float(kappa), 2),
                round(float(alpha95), 1),
                round(float(eigenvalues[0]), 6),
                round(float(eigenvalues[1]), 6),
                round(float(eigenvalues[2]), 6),
                round(float(az_eigen1), 1),
                round(float(inc_eigen1), 1),
                round(float(az_eigen2), 1),
                round(float(inc_eigen2), 1),
                round(float(az_eigen3), 1),
                round(float(inc_eigen3), 1),
            ])
            mean_bedding_provider.addFeature(mean_point)
                        
            wc_names.append(polygon_name)
            wc_eigenvalues.append([eigenvalues[0], eigenvalues[1], eigenvalues[2]])
            wc_n.append(n)
    
    
    mean_bedding_layer.commitChanges()

    # Add unused points to the 'Unused_points' layer
    for point in bedding_layer.getFeatures():
        if point.id() not in used_points:
            unused_feature = QgsFeature()
            unused_feature.setGeometry(point.geometry())
            unused_feature.setAttributes(point.attributes())
            unused_provider.addFeature(unused_feature)
            
    # Ensure that points that are alone in polygons are also saved
    for area_feature in areas_layer.getFeatures():
        geom = area_feature.geometry()
        points_inside = [point_feature for point_feature in all_points_inside if geom.contains(point_feature.geometry())]

        # If there is only one point in the polygon
        if len(points_inside) == 1:
            point_feature = points_inside[0]
            unused_feature = QgsFeature()
            unused_feature.setGeometry(point_feature.geometry())
            unused_feature.setAttributes([
                point_feature.geometry().asPoint().x(),
                point_feature.geometry().asPoint().y()
            ])
            unused_provider.addFeature(unused_feature)  # Add to the unused points layer
    
    

     # Add layers to the project
    QgsProject.instance().addMapLayer(mean_bedding_layer)
    QgsProject.instance().addMapLayer(unused_points_layer)
        
    # After processing, save the table to CSV
    project_path = QgsProject.instance().homePath()  # Get the project directory
    csv_path = os.path.join(project_path, 'Mean_Bedding.csv')  # Define the path for the CSV


    # Save the layer in CSV format
    error = QgsVectorFileWriter.writeAsVectorFormat(
        mean_bedding_layer, 
        csv_path, 
        'UTF-8', 
        mean_bedding_layer.crs(), 
        'CSV'
    )

    if error[0] == QgsVectorFileWriter.NoError:
        print(f'Archivo guardado correctamente en {csv_path}')
    else:
        print('Error al guardar el archivo:', error[1])
        
    # Drawing the Woodcock plot
    woodcock_plot(wc_names, wc_eigenvalues, wc_n)
    

mean_bedding_main()
