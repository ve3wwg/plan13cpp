The BASIC program and info was obtained from the following site:

    https://www.amsat.org/articles/g3ruh/111.html
    (C)1990 J.R. Miller G3RUH

Plan10 and then Plan13 was developed by James Miller G3RUH, from
OSCAR-10 position Calculation Program, Oscar News, 1983 Dec, No.45
p.15-18.

The C++ code was downloaded from:

    https://github.com/BackupGGCode/qrptracker.git

I couldn't determine who transcribed it to Arduino (C++). I have
modified it to be C++/Linux friendly, removing the Arduino
references.

Modify main.cpp to reflect your own location and recent keplerian
elements. 

For Linux/*BSD:

    $ make
    $ ./main 

To test:

    $ wget -q http://api.open-notify.org/iss-now.json -O - && echo && ./main
    {"timestamp": 1567815574, "iss_position": {"longitude": "4.9295", "latitude": "39.1889"}, "message": "success"}
      Date / Time UTC    Sat   Azimuth   Elevation   Latitude   Longitude   RR
    ------------------- ----- ---------- ---------- ---------- ---------- ------
    2019-09-07 00:19:34 ISS    62.032916 -27.088588   39.23168    4.94131  1.914
    2019-09-07 00:19:34 AO-7   68.543398 -55.327083   -4.81282   45.17617  1.713
    2019-09-07 00:19:34 UO-11 355.884043 -47.113599   37.75230  106.15994  4.657

and compare the ISS Latitude and Longitude from the JSON to the printed
output.

There was no license file found at the above githup repo.

My own changes are public domain.
