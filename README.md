## Numerical Temperature prediction
===================================

Creators: Tayfun Kraderi, Falola Yusuf

This software is a C++ program created to numerically predict tempeartures in Heathrow, London. We ensured at least 50 data points were used for each parameter in this study. The average temperature, that is, the average value of minimum and maximum temperature, was considered. Parameters considered to predict temperature include; year, CO2 emission, sunspot number, and population growth.


## Numerical methods adopted:
- Single variable linear regression


- Multivariable linear regression 

## Usage
- GitHub page should be cloned, and all files saved in the folder of the new project created.

- The user should open the Temperature_Forecast.vcxproj file. It is assumed the user has MSVC 2017 installed. Then folder directory should be included by right clicking on Temperature_Forecast.vcxproj >> properties >> C/C++ >> General >> Additional Include Directories in MSVC 2017.

- Once the program is run, the user would be prompted to enter the year for which he/she wants the temperature to be predicted, then the month of prediction, (1 for January and 12 for December and so on).

- The temperature prediction for the two models used would show up for the specified month and year.

- This result would also be automatically saved in the prediction.txt file.


## Data
- The data.txt file contains the data (attributes) used for temperature prediction with columns as shown below;

   Year | Population| Sunspot number| CO2 emission 
   
- The temperature.txt files contains historical temperature data with columns as follows;

   Year | Month| Min Temperature (deg C)| Max Temperature (deg C)| af (days) | Rain (mm) | Sun (hrs)

## Data sources:
1. https://data.london.gov.uk/dataset/office-national-statistics-ons-population-estimates-borough?resource=20dc1341-e74a-4e20-b1ff-a01c45e9fa10

2. http://www.sidc.be/silso/datafiles-old

3. https://data.worldbank.org/indicator/EN.ATM.CO2E.PC

