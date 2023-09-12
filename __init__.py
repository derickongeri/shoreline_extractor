# -*- coding: utf-8 -*-
"""
/***************************************************************************
 AutomaticShorelineExtraction
                                 A QGIS plugin
 This plugin enables the user to automatically extract shorelines and compute shoreline change rates.
 Generated by Plugin Builder: http://g-sherman.github.io/Qgis-Plugin-Builder/
                             -------------------
        begin                : 2023-09-07
        copyright            : (C) 2023 by LocateIT
        email                : cogeos@locateit.co.ke
        git sha              : $Format:%H$
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
 This script initializes the plugin, making it known to QGIS.
"""


# noinspection PyPep8Naming
def classFactory(iface):  # pylint: disable=invalid-name
    """Load AutomaticShorelineExtraction class from file AutomaticShorelineExtraction.

    :param iface: A QGIS interface instance.
    :type iface: QgsInterface
    """
    #
    from .shoreline_extractor import AutomaticShorelineExtraction
    return AutomaticShorelineExtraction(iface)
