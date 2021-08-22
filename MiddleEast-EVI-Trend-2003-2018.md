# EVI Change Trends in the Middle East from 2003 to 2018:
I created this map to tell the region stories by capture and emphasize the long-term slop of [EVI](https://en.wikipedia.org/wiki/Enhanced_vegetation_index)
using the [harmonic regression](https://docs.google.com/document/d/1mNIRB90jwLuASO1JYas1kuOXCLbOoy1Z4NlV1qIXM10/edit) model calculated via 
[Google Earth Engine](https://developers.google.com/earth-engine/guides/playground). 
The used data source is the [MOD13Q1.006](https://developers.google.com/earth-engine/datasets/catalog/MODIS_006_MOD13Q1) Terra Vegetation Indices 
16-Day Global 250m. Also, you can test it online [here](https://code.earthengine.google.com/9508fccb9c007c50ddb86c91ffe1f2c7).

Here is a picture that worth a thousand words:
![evi trend stories](https://user-images.githubusercontent.com/11270404/130354551-3d1d7e6c-d251-468f-93e6-048cd1184661.jpg)

> Models Uncertainty: Essentially, all models are wrong, but some are useful ([George E. P. Box](https://en.wikipedia.org/wiki/George_E._P._Box)).

* 1<sup>st</sup> Story: [Urban Expansion vs. Land Reclaiming](#1st-story-urban-expansion-vs-land-reclaiming)
* 2<sup>nd</sup> Story: [Pivot-Irrigated Fields in the Desert](#2nd-story-pivot-irrigated-fields-in-the-desert)
* 3<sup>rd</sup> Story: [Drought and Water Level](#3rd-story-drought-and-water-level)
* 4<sup>th</sup> Story: [Crises Impact on Agriculture](#4th-story-crises-impact-on-agriculture)

Find below the script that I developed to generate this map:

```javascript
// Library that includes all the utility functions for ease of use
var Lib = require("users/khaledalshamaa/icarda:Modules/MOD13Q1TrendModule.js");

// Import MOD13Q1.006 Terra Vegetation Indices 16-Day Global 250m
// https://code.earthengine.google.com/dataset/MODIS/006/MOD13Q1
var MOD13Q1 = ee.ImageCollection("MODIS/006/MOD13Q1");

// Import LSIB 2017: Large Scale International Boundary Polygons, Simplified
// https://developers.google.com/earth-engine/datasets/catalog/USDOS_LSIB_SIMPLE_2017
var LSIB = ee.FeatureCollection("USDOS/LSIB_SIMPLE/2017");
var ROI  = LSIB.filter(ee.Filter.inList("country_na", Lib.MiddleEast));

// Setup the base map
Map.setOptions("Gray", {"Gray": Lib.grayMapStyle});
Map.addLayer(ROI, {}, "Middle East");
Map.centerObject(ROI, 5);

// Add the legend
var palette = ["529400","74A901","FFFFFF","FAB505","DF923D"];
var names = ["+0.2","+0.1"," 0.0","-0.1","-0.2"];
 
Map.add(Lib.legend(palette, names, "EVI Changes", "bottom-left"));

// Study period start-end (yyyy-mm-dd)
var start = "2003-01-01";
var end   = "2018-12-31";

var period = ee.Date(end).difference(ee.Date(start), "year");

// Time and Spatial Filtering
var MOD13Q1ME = MOD13Q1.filterDate(start, end).filterBounds(ROI);

// Apply addNDVI on all images in the MODIS collection 
MOD13Q1ME = MOD13Q1ME.map(Lib.maskModisVI).map(Lib.addParams);

// Estimate seasonality with a harmonic model
var harmonicModel = MOD13Q1ME.map(Lib.addHarmonicTerms);

// Name of the dependent variable
var dependent = ee.String("EVI_rescaled");

// Define the independent variables for the harmonic regression
var independents = Lib.getHarmonicIndependents();

// The output of the regression reduction is a 4x1 array image.
var trend = harmonicModel.select(independents.add(dependent))
                         .reduce(ee.Reducer.linearRegression(independents.length(), 1));

// Turn the array image into a multi-band image of coefficients.
var trendCoefficients = trend.select("coefficients")
                             .arrayProject([0])
                             .arrayFlatten([independents]);

// Get the long term slop 
var longTermChanges = trendCoefficients.select("t").multiply(period).clip(ROI);

Map.addLayer(longTermChanges, 
            {
              min: Number(names[names.length - 1]), 
              max: Number(names[0]), 
              palette: palette.reverse().join(",")
            }, 
            "Long Term EVI Changes");

// Compute the fitted (predicted) values
var fittedHarmonic = function(image) {
  var predicted = image.select(independents)
                       .multiply(trendCoefficients)
                       .reduce("sum")
                       .rename("fitted");

  return image.addBands(predicted);
};

harmonicModel = harmonicModel.map(fittedHarmonic);

// Plot the fitted model and the original data at the Location Loc
function getChart(coords){
  var point = ee.Geometry.Point(coords.lon, coords.lat).buffer(250);
  var fittedVsReal = ui.Chart.image.series(harmonicModel.select(["fitted","EVI_rescaled"]), 
                                           point, ee.Reducer.mean(), 1000);

  fittedVsReal.setSeriesNames(["EVI", "fitted"]);
  fittedVsReal.style().set({width: "400px", position: "bottom-right"});

  Map.widgets().set(1, fittedVsReal);
}

// Add interactive to show the time series by select point on a map
Map.onClick(getChart);
Map.style().set("cursor", "crosshair");
```

**_Author:_** Khaled Al-Shamaa <k.el-shamaa@cgiar.org>

**_License:_** This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by 
the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

**_Disclaimer:_** This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY 
or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

---

### Normalized Difference Vegetation Index (NDVI):
[NDVI](https://earthobservatory.nasa.gov/Features/MeasuringVegetation/measuring_vegetation_2.php) quantifies vegetation by measuring the difference between 
near-infrared (which vegetation strongly reflects) and red light (which vegetation absorbs):

`NDVI = (NIR - RED) / (NIR + RED)`

NDVI always ranges from -1 to +1 (water NDVI value is negative). Healthy vegetation absorbs the most visible light that reaches it and reflects a great deal 
of near-infrared light. It has evolved to do this because absorbing more infrared light would overheat the plants.

### Enhanced Vegetation Index (EVI):
In areas of the dense canopy where the leaf area index (LAI) is high, the NDVI values can be improved by leveraging information in the blue wavelength. 
Information in this portion of the spectrum can help correct for soil background signals and atmospheric influences.

`EVI = 2.5 x (NIR - RED) / (NIR + 6 x RED - 7.5 x BLUE + 1)`

The range of values for the EVI is -1 to 1, where healthy vegetation generally falls between values of 0.20 to 0.80.

---

## 1st Story: Urban Expansion vs. Land Reclaiming
![story1](https://user-images.githubusercontent.com/11270404/130354565-315b6f39-af61-4c08-9dae-ab7838ff322a.jpg)

> Reference point(s): urban sample 31.57E, 30.40N & reclamation sample 31.57E, 30.35N

## 2nd Story: Pivot-Irrigated Fields in the Desert
![story2](https://user-images.githubusercontent.com/11270404/130354592-470015a2-9bc6-40e0-893b-72186da19805.jpg)

> Reference point(s): 38.3E, 30.0N

## 3rd Story: Drought and Water Level
![story3](https://user-images.githubusercontent.com/11270404/130354599-8360544a-fe6a-4cf6-95e6-61bb3cdafb8f.jpg)

> Reference point(s): 43.41E, 34.06N

## 4th Story: Crises Impact on Agriculture
![story4](https://user-images.githubusercontent.com/11270404/130354609-0e349416-f897-4b65-bf26-e28c2e04b5dd.jpg)

> Reference point(s): Turkey sample 40.33E, 37.00N & Syria sample 40.48E, 36.92N

