import sami2py

day = 208
year = 2012

sami2py.run_model(day=day,year=year,hrmax=24.5,fmtout=True,tag='test')

S = sami2py.Sami2Model(tag='test',lon=0,year=year,day=day)
