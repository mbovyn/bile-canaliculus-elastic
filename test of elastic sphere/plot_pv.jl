
using CSV, DataFrames, PlotlyJS

#set working directory to directory of script
cd(@__DIR__)

# take in data written out by surface evolver
#linear_elastic, original try (two refinements)
df = CSV.File("sphere_pV.txt") |> DataFrame
#linear_elastic, two refinements
df2 = CSV.File("sphere_pV_r2.txt") |> DataFrame
#linear_elastic, three refinements
df3 = CSV.File("sphere_pV_r3.txt") |> DataFrame
#linear_elastic, four refinements
df4 = CSV.File("sphere_pV_r4.txt") |> DataFrame
#general_linear_elastic, two refinements
dfgen = CSV.File("sphere_pV_gen.txt") |> DataFrame

## Make a basic plot of the data

layout = Layout(
    title="Surface evolver sphere",
    xaxis_title="Pressure",
    yaxis_title="Volume"
)

plot(df,x=:p,y=:v,layout)

## Theory curves

#pressures
thp=range(0,stop=1,length=100)
#material parameters
e=1
h0=1
nu=.5
#get initial radius from data
r0=cbrt(3/(4*pi)*df.v[1])
#volume from Knoche and Kierfeld paper (KK theory)
kkv=4/3*pi*(e*h0./(thp*(1-nu))-sqrt.((e*h0./(thp*(1-nu))).^2-(2*e*h0)./(thp*(1-nu))*r0)).^3
#volume from the equation Pierre wrote down
phv=4/3*pi*(thp.*r0^2 ./(4*e*h0).+ r0*sqrt.(1 .+((thp.*r0)./(4*e*h0)).^2)).^3
#volume from the equation I got in Mathematica
mbv=4/3*pi*( -(nu-1)*r0/(2*e*h0) * ( thp*r0 .+ sqrt.(thp.^2*r0^2 .+ (4*(e*h0)^2/(nu-1)^2)) ) ).^3

## Make plots

#first try at theory plot
#plot(scatter(x=thp,y=kkv))

layout2 = Layout(
    title="Elastic sphere",
    xaxis_title="Pressure",
    yaxis_title="Volume",
)

#Different curves to plot
#linear_elastic, original try (2 refinements) (several ways)
#s=scatter(df,x=:p,y=:v,name="Surface Evolver")
#s=scatter(x=df.p.-df.p[1],y=df.v)
s2=scatter(df2,x=:p,y=:v,name="Surface Evolver")
#linear_elastic, two refinements
#s2=scatter(df2,x=:p,y=:v,name="SE, r2")
#linear_elastic, three refinements
s3=scatter(df3,x=:p,y=:v,name="SE, r3")
#linear_elastic, four refinements
s4=scatter(df4,x=:p,y=:v,name="SE, r4")
#general_linear_elastic, two refinements
s5=scatter(dfgen,x=:p,y=:v,name="general_linear_elastic")
#linear_elastic, original try (2 refinements)
s6=scatter(df,x=:p,y=:v,name="linear_elastic")

#KK theory and original try
#p=plot([s,scatter(x=thp,y=kkv,name="Theory")],layout2)
#KK theory and different numbers of refinements
p2=plot([s2,s3,s4,scatter(x=thp,y=kkv,name="e=λ-1")],layout2)
#KK and "non-linear" theory with linear and general
p=plot([s6,s5,scatter(x=thp,y=mbv,name="e=(λ²-1)/2"),scatter(x=thp,y=kkv,name="e=λ-1")],layout2)

## save figure
savefig(p,"elastic sphere result linearizations.pdf");
savefig(p2,"elastic sphere result linearizations.pdf");