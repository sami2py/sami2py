import sami2py

day = 210
year = 2012
lon = 0

sami2py.run_model(day=day,year=year,lon=lon,hrmax=24.5,
                  fmtout=True,tag='test',test=True,o_scale=0.8)

S = sami2py.model(tag='test',lon=lon,year=year,day=day)
