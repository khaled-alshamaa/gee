# GEE Script Calculates EVI Seasonal Peak:
I created this beautiful map that shows the seasonal [EVI](https://en.wikipedia.org/wiki/Enhanced_vegetation_index) peak and magnitude during 
specific period using [Google Earth Engine](https://developers.google.com/earth-engine/guides/playground). 
It calculates the [harmonic regression](https://docs.google.com/document/d/1mNIRB90jwLuASO1JYas1kuOXCLbOoy1Z4NlV1qIXM10/edit) of the EVI time series, 
then translates the phase and amplitude of the estimated model into [HSV](https://en.wikipedia.org/wiki/HSL_and_HSV) color space. 
Please note that presented amplitude (i.e., color saturation) refers to the seasonal component not the absolute EVI value. 

The used data source is the [MOD13Q1.006](https://developers.google.com/earth-engine/datasets/catalog/MODIS_006_MOD13Q1) Terra Vegetation Indices 
16-Day Global 250m. Here is the script that I developed to generate this [map](https://code.earthengine.google.com/46d407e057c09eb6edb8bf2f57c67fe0).

![Seasonal EVI Peak](https://user-images.githubusercontent.com/11270404/130336574-bd096dc3-85ad-48ca-9523-6fe0e7761d3e.png)

## Required Utility Functions:
This function receives a MOD13Q1 image and returns same image with an extra band for time in fractional years since the epoch (1970-01-01) name it as "t" 
and "intercept" band of constant value 1 will be used for regression analysis.

```javascript
var addParams = function(img){
  var ndvi = img.select("NDVI").multiply(0.0001).rename("NDVI_rescaled");
  var evi  = img.select("EVI").multiply(0.0001).rename("EVI_rescaled");
  
  // compute time in fractional years since the epoch
  var date  = ee.Date(img.get("system:time_start"));
  var years = date.difference(ee.Date("1970-01-01"), "year");
  
  // generate time and constant (i.e. intercept) bands for regression model
  var t = ee.Image(years).rename("t").float();
  var a = ee.Image.constant(1).rename("intercept");
  
  return img.addBands(ndvi).addBands(evi).addBands(t).addBands(a);
};
```

Mask upon VI quality indicator(s):
```javascript
var maskModisVI = function(image){
  var qa   = image.select("DetailedQA");
  var mask = qa.bitwiseAnd(2).eq(0);
  
  return image.updateMask(mask);
};
```

Estimate seasonality with a harmonic model:
> pt = β0 + β1*t + β2*cos(2πωt) + β3*sin(2πωt) + et
 
where w refers to how many cycles per unit time (i.e., year).
 
```javascript
var w = 2;

var addHarmonicTerms = function(image){
  for(var i = 1; i <= w; i++){
    // calculate the term: 2πωt
    var timeRadians = image.select("t").multiply(2 * i * Math.PI);
    
    // calculate the terms cos(2πωt) and sin(2πωt)
    var cos = timeRadians.cos().rename("cos"+i);
    var sin = timeRadians.sin().rename("sin"+i);
    
    image = image.addBands(cos).addBands(sin);
  }
  
  return image;
};
```
 
Define the independent variables for the harmonic regression:
```javascript
var getHarmonicIndependents = function(){
  var independents = ee.List(["intercept", "t"]);
  
  for(var i = 1; i <= w; i++){
    independents = ee.List(independents).add("cos"+i);
    independents = ee.List(independents).add("sin"+i);
  }
  
  return independents;
};
```
 
Create a customized legend block:
```javascript
var legend = function(palette, names, lbl, pos) {
  // set position of panel
  var legend = ui.Panel({style:{position:pos, padding:"8px 15px"}});
   
  // create legend title
  var legendTitle = ui.Label({value:lbl, style:{margin:"0 0 4px 0", padding:"0"}});
   
  // add the title to the panel
  legend.add(legendTitle);
   
  // creates and styles 1 row of the legend.
  var makeRow = function(color, name) {
    // create the label that is actually the colored box.
    var colorBox = ui.Label({style:{backgroundColor:"#"+color, border:"1px solid black", padding:"8px", margin:"0 0 4px 0"}});
  
    // create the label filled with the description text.
    var description = ui.Label({value:name, style:{margin:"0 0 4px 6px"}});
  
    // return the panel
    return ui.Panel({widgets:[colorBox, description], layout:ui.Panel.Layout.Flow("horizontal")});
  };

  // add color and and names
  for (var i = 0; i < palette.length; i++) legend.add(makeRow(palette[i], names[i]));
  
  return legend;
};
```
 
## Main Script Code:
 
```javascript
var start = "2018-01-01";
var end   = "2020-12-31";

// import MOD13Q1.006 Terra Vegetation Indices 16-Day Global 250m
var modis = ee.ImageCollection("MODIS/006/MOD13Q1");

// temporal filtering
var modis = modis.filterDate(start, end);

// mask upon QA and add required regression parameters
modis = modis.map(maskModisVI).map(addParams);

// calculate harmonic regression terms
var harmonicModel = modis.map(addHarmonicTerms);

// name of the dependent variable
var dependent = ee.String("EVI_rescaled");

// define the independent variables for the harmonic regression
var independents = getHarmonicIndependents();

// the output of the regression reduction is a 4x1 array image.
var trend = harmonicModel.select(independents.add(dependent)).reduce(ee.Reducer.linearRegression(independents.length(), 1));

// turn the array image into a multi-band image of coefficients.
var trendCoefficients = trend.select("coefficients").arrayProject([0]).arrayFlatten([independents]);

// map phase and amplitude of the estimated harmonic model
// compute phase and amplitude 
// atan2: calculates the angle formed by the 2D vector [x, y]
// hypot: calculates the magnitude of the 2D vector [x, y]
var phase     = trendCoefficients.select("cos1").atan2(trendCoefficients.select("sin1"));
var amplitude = trendCoefficients.select("cos1").hypot(trendCoefficients.select("sin1"));

// use the HSV to RGB transform to display phase and amplitude
var rgb = phase.unitScale(-Math.PI, Math.PI).addBands(amplitude.multiply(2.5)).addBands(ee.Image(1)).hsvToRgb();

Map.addLayer(rgb, {}, "Phase (hue), Amplitude (saturation)");

// palette with the colors
var palette = ["00BBFF","0044FF","4400FF","CC00FF","FF00BB","FF0033","FF4400","FFCC00","BBFF00","33FF00","00FF44","00FFCC"];
 
// name of the legend
var names = ["Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"];

Map.add(legend(palette, names, "EVI Peak", "bottom-left"));
```

**_Author:_** Khaled Al-Shamaa <k.el-shamaa@cgiar.org>

**_License:_** This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

**_Disclaimer:_** This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
