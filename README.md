# Characteristics and trends of Australian climate-change policies

<img align="right" src="www/climchangepolicy.jpg" alt="climate change" width="180" style="margin-top: 20px">

<br>
Prof <a href="https://globalecologyflinders.com/people/#DIRECTOR">Corey J. A. Bradshaw</a> <br>
<a href="http://globalecologyflinders.com" target="_blank">Global Ecology</a> | <em><a href="https://globalecologyflinders.com/partuyarta-ngadluku-wardli-kuu/" target="_blank">Partuyarta Ngadluku Wardli Kuu</a></em>, <a href="http://flinders.edu.au" target="_blank">Flinders University</a>, Adelaide, Australia <br>
April 2026 <br>
<a href=mailto:corey.bradshaw@flinders.edu.au>e-mail</a> <br>
<br>
Maddy King
<a href="http://flinders.edu.au" target="_blank">Flinders University</a>, Adelaide, Australia <br>
<a href=mailto:maddy.king@flinders.edu.au>e-mail</a> <br>
<br>
Dr <a href="https://www.flinders.edu.au/people/altaf.virani">Altaf Virani</a><br>
<a href="[http://globalecologyflinders.com](https://www.flinders.edu.au/college-business-creative-arts-law-social-sciences)" target="_blank">College of Business, Government and Law</a><br>
<a href="http://flinders.edu.au" target="_blank">Flinders University</a>, Adelaide, Australia <br>
<a href=mailto:altaf.virani@flinders.edu.au>e-mail</a> <br>
<br>

## <a href="https://github.com/cjabradshaw/ClimateChangeAdaptationPolicies/tree/main/scripts">Scripts</a>
- <code>CAP.R</code> (main code)
- <code>SILOncDownload.txt</code>: Terminal commands to download annual .nc layers for Australia from <a href="https://www.longpaddock.qld.gov.au/silo/">SILO</a>

## Data
- <em>CAPdat.csv</em>: main policy database
- <em>earningMedSA4</em>: median individual earnings per <a href="https://www.abs.gov.au/statistics/standards/australian-statistical-geography-standard-asgs-edition-3/jul2021-jun2026/main-structure-and-greater-capital-city-statistical-areas/statistical-area-level-4">Statistical Area Level 4</a> (SA4) (<a href="https://www.abs.gov.au/statistics/labour/earnings-and-working-conditions/personal-income-australia/2021-22#data-downloads">Australian Bureau of Statistics</a>)
- annual NetCDF (.nc) climate layers for Australia downloaded by the user (see <code>SILOncDownload.txt</code> in the '<a href="https://github.com/cjabradshaw/ClimateChangeAdaptationPolicies/tree/main/scripts">scripts</a>' subdirectory)
- <em>popSA4.csv</em>: resident population estimates by SA4 code from 2001–2025 (<a href="https://dataexplorer.abs.gov.au/vis?tm=ABS_ANNUAL_ERP_ASGS2021&pg=0&snb=1&df%5Bds%5D=PEOPLE_TOPICS&df%5Bid%5D=ABS_ANNUAL_ERP_ASGS2021&df%5Bag%5D=ABS&df%5Bvs%5D=1.2.0&dq=.GCCSA..A&pd=2015,&to%5BTIME_PERIOD%5D=false">Australian Bureau of Statistics</a>)

## R libraries
<code>data.table</code>, <code>dplyr</code>, <code>ggplot2</code>, <code>ggpubr</code>, <code>ggrepel</code>, <code>lubridate</code>, <code>ozmaps</code>, <code>purrr</code>, <code>rnaturalearth</code>, <code>scales</code>, <code>sf</code>, <code>terra</code>, <code>tidyr</code>, <code>units</code>, <code>viridis</code>
<br>
<br>

<p><a href="https://www.flinders.edu.au"><img align="bottom-left" src="www/Flinders_University_Logo_Stacked_RGB_Master.jpg" alt="Flinders University" width="80" style="margin-top: 20px"></a> &nbsp; <a href="https://globalecologyflinders.com"><img align="bottom-left" src="www/GEL Logo Kaurna New Transp.png" alt="GEL" width="170" style="margin-top: 20px"></a></p>

