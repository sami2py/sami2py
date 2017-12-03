import sami2py

day = 209
year = 2012

sami2py.run_model(day=day,year=year,hrmax=24.5,fmtout=False,tag='test')

S = sami2py.model(tag='test',lon=0,year=year,day=day)
