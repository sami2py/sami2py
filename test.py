import sami2py

day = 211
year = 2012
lon = 0

#sami2py.run_model(day=day,year=year,lon=lon,hrmax=25.0,
#                  fmtout=True,tag='test')

S = sami2py.model(tag='test',lon=lon,year=year,day=day)
