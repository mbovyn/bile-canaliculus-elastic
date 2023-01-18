using CSV, DataFrames, PlotlyJS, DataInterpolations

cd(@__DIR__)
nfiles=2

scl=3.6

## get p-v curves

global sss=String[]
for i=0:nfiles-1
    push!(sss, string("p-v/","v-",i,".txt"))
end

dfs = DataFrame.(CSV.File.(sss))

vcrit=[5,4.85096]

## diam v volume


plot([scatter(dfs[1],y=:maxx,x=:v),
    scatter(dfs[2],y=:maxx,x=:v)]
)

# x0interp=LinearInterpolation(dfs[1].maxx,dfs[1].v)
# x0interp(vcrit[1])

# x1interp=LinearInterpolation(dfs[2].maxx,dfs[2].v)
# x1interp(vcrit[2])

icrit=3
xcrit=[dfs[1].maxx[icrit],dfs[2].maxx[icrit]]


## Normalized version

plot([scatter(dfs[1],y=dfs[1].maxx/dfs[1].maxx[icrit],x=dfs[1].v/dfs[1].v[icrit]),
    scatter(dfs[2],y=dfs[2].maxx/dfs[2].maxx[icrit],x=dfs[2].v/dfs[2].v[icrit])]
)

x0interp=LinearInterpolation(dfs[1].maxx/dfs[1].maxx[icrit],dfs[1].v/dfs[1].v[icrit])
x0interp(1.76)

x1interp=LinearInterpolation(dfs[2].maxx/dfs[2].maxx[icrit],dfs[2].v/dfs[2].v[icrit])
x1interp(1.76)

v0interp=LinearInterpolation(dfs[1].v/dfs[1].v[icrit],dfs[1].maxx/dfs[1].maxx[icrit])
v1interp=LinearInterpolation(dfs[2].v/dfs[2].v[icrit],dfs[2].maxx/dfs[2].maxx[icrit])
[v0interp(1.2),v1interp(1.2)]


## Version with pressure as a function of volume

pvp=plot(Layout(
    xaxis_title="Volume (μm³)",
    yaxis_title="Pressure/Cortex Modulus",
    width=400,
    height=400,
    legend=attr(x=0.025, y=.975),
    xaxis=attr(tickmode="auto",nticks=10),
    xaxis_range=[170,400]
))
#pvp

fac=v1interp(1.2)

vpinterp1=LinearInterpolation(dfs[2].p,dfs[2].v)
p1f=vpinterp1(vcrit[2]*fac)
p10=vpinterp1(vcrit[2])
dp1=p1f-p10

vpinterp0=LinearInterpolation(dfs[1].p,dfs[1].v)
p0f=vpinterp0(vcrit[1]*fac)
p00=vpinterp0(vcrit[1])
dp0=p0f-p00

pmore=dp1/dp0
print(pmore)

rightlim=10

addtraces!(pvp,
scatter(
    x=[vcrit[2],vcrit[2]*fac,rightlim,rightlim]*scl^3,
    y=[p10,p1f,p1f,p10],
    fill="toself",mode="lines",fillcolor="Red",showlegend=false,opacity=.2,line_color=:transparent),
scatter(
    x=[vcrit[1],vcrit[1]*fac,rightlim,rightlim]*scl^3,
    y=[p00,p0f,p0f,p00],
    fill="toself",mode="lines",fillcolor="RoyalBlue",showlegend=false,opacity=.2,line_color=:transparent),
)

addtraces!(pvp,
    scatter(dfs[2],x=dfs[2].v*scl^3,y=:p,name="With Bulkheads",line_color="Red"),
    scatter(dfs[1],x=dfs[1].v*scl^3,y=:p,name="No Bulkheads",line_color="RoyalBlue"),
)
pvp

## version with pressure as a function of diameter

pw=400
ph=350

pvd=plot(Layout(
    xaxis_title="Diameter (μm)",
    yaxis_title="Pressure/Cortex Modulus",
    width=pw,
    height=ph,
    legend=attr(x=.05, y=.85),
    xaxis=attr(tickmode="auto",nticks=7),
    xaxis_range=[3.25,5],
    yaxis_range=[0,2.75]
))
#pvp

fac=1.21

xpinterp1=LinearInterpolation(dfs[2].p,dfs[2].maxx)
p1f=xpinterp1(xcrit[2]*fac)
p10=xpinterp1(xcrit[2])
dp1=p1f-p10

xpinterp0=LinearInterpolation(dfs[1].p,dfs[1].maxx)
p0f=xpinterp0(xcrit[1]*fac)
p00=xpinterp0(xcrit[1])
dp0=p0f-p00

pmore=dp1/dp0
print(pmore)

rightlim=10

addtraces!(pvd,
scatter(
    x=[xcrit[2],xcrit[2]*fac,rightlim,rightlim]*scl,
    y=[p10,p1f,p1f,p10],
    fill="toself",mode="lines",fillcolor="Green",showlegend=false,opacity=.2,line_color=:transparent),
scatter(
    x=[xcrit[1],xcrit[1]*fac,rightlim,rightlim]*scl,
    y=[p00,p0f,p0f,p00],
    fill="toself",mode="lines",fillcolor="Black",showlegend=false,opacity=.2,line_color=:transparent),
)

addtraces!(pvd,
    scatter(dfs[2],x=dfs[2].maxx*scl,y=:p,name="With Bulkheads",line_color="Green"),
    scatter(dfs[1],x=dfs[1].maxx*scl,y=:p,name="No Bulkheads",line_color="Black"),
)

addtraces!(pvd,
    scatter(x=[1,1]*xcrit[2]*scl,y=[dfs[2].p[icrit]+.1,2.6],text=["","D*"],mode="lines+text",textposition="left",textfont=attr(color="#6a9c7a"),line_color="#6a9c7a",line=attr(dash="dot"),showlegend=false),
    scatter(x=[1,1].*xcrit[2]*fac*scl,y=[p1f+.1,2.6],text=["","Dbleb "],mode="lines+text",textposition="left",textfont=attr(color="#6a9c7a"),line_color="#6a9c7a",line=attr(dash="dot"),showlegend=false),
    scatter(x=[1,1]*xcrit[1]*scl,y=[dfs[1].p[icrit]-.1,.1],text=[""," D*"],mode="lines+text",textposition="right",textfont=attr(color="Black"),line_color="Grey",line=attr(dash="dot"),showlegend=false),
    scatter(x=[1,1].*xcrit[1]*fac*scl,y=[p0f,.1],text=[""," Dbleb"],mode="lines+text",textposition="right",textfont=attr(color="Black"),line_color="Grey",line=attr(dash="dot"),showlegend=false),
)

addtraces!(pvd,
    scatter(x=[1,1]*xcrit[2]*scl,y=[dfs[2].p[icrit]+.1,2.6],text=["","D*"],mode="lines+text",textposition="left",textfont=attr(color="#6a9c7a"),line_color="#6a9c7a",line=attr(dash="dot"),showlegend=false),
    scatter(x=[1,1].*xcrit[2]*fac*scl,y=[p1f+.1,2.6],text=["","Dbleb "],mode="lines+text",textposition="left",textfont=attr(color="#6a9c7a"),line_color="#6a9c7a",line=attr(dash="dot"),showlegend=false),
    scatter(x=[1,1]*xcrit[1]*scl,y=[dfs[1].p[icrit]-.1,.1],text=[""," D*"],mode="lines+text",textposition="right",textfont=attr(color="Black"),line_color="Grey",line=attr(dash="dot"),showlegend=false),
    scatter(x=[1,1].*xcrit[1]*fac*scl,y=[p0f,.1],text=[""," Dbleb"],mode="lines+text",textposition="right",textfont=attr(color="Black"),line_color="Grey",line=attr(dash="dot"),showlegend=false),
)
pvd

##
#savefig(pvd,"pV_and_deltap_new.pdf",width=pw,height=ph)

##

dfs[2].p./xpinterp0.(dfs[2].maxx)

## New version focusing on pressure difference at constant diameter

pscl=400/300
pw=ceil(Int64,325*pscl)
ph=ceil(Int64,250*pscl)

pvd2=plot(Layout(
    xaxis_title="Diameter (μm)",
    yaxis_title="Pressure/Cortex Modulus",
    width=pw,
    height=ph,
    #legend=attr(x=.01, y=.89),
    xaxis=attr(tickmode="auto",nticks=7),
    xaxis_range=[3.25,5],
    yaxis_range=[0,2.75].*.2,
    font=attr(size=8*pscl,family="Arial")
))

addtraces!(pvd2,
scatter(
    x=[xcrit[2],xcrit[2],rightlim,rightlim]*scl,
    y=[p10,p00,p00,p10].*.2,
    fill="toself",mode="lines",fillcolor="Green",showlegend=false,opacity=.2,line_color=:transparent),
)

addtraces!(pvd2,
scatter(
    x=[5],
    y=[(p00+p10)/2].*.2,
    text=["Δp"],textfont=attr(family="Arial"),mode="text",textposition="left",showlegend=false),
)

addtraces!(pvd2,
    scatter(dfs[2],x=dfs[2].maxx*scl,y=dfs[2].p.*.2,name="With Bulkheads",line_color="Green"),
    scatter(dfs[1],x=dfs[1].maxx*scl,y=dfs[1].p.*.2,name="No Bulkheads",line_color="Black"),
)

addtraces!(pvd2,
    scatter(x=[1,1]*xcrit[2]*scl,y=[dfs[2].p[icrit]+.1,2.6].*.2,text=["","D*"],mode="lines+text",textposition="left",textfont=attr(color="#6a9c7a"),line_color="#6a9c7a",line=attr(dash="dot"),showlegend=false),
    scatter(x=[1,1].*xcrit[2]*fac*scl,y=[p1f+.1,2.6].*.2,text=["","Dbleb "],mode="lines+text",textposition="left",textfont=attr(color="#6a9c7a"),line_color="#6a9c7a",line=attr(dash="dot"),showlegend=false),
    scatter(x=[1,1]*xcrit[1]*scl,y=[dfs[1].p[icrit]-.1,.1].*.2,text=[""," D*"],mode="lines+text",textposition="right",textfont=attr(color="Grey"),line_color="Grey",line=attr(dash="dot"),showlegend=false),
    scatter(x=[1,1].*xcrit[1]*fac*scl,y=[p0f,.1].*.2,text=[""," Dbleb"],mode="lines+text",textposition="right",textfont=attr(color="Grey"),line_color="Grey",line=attr(dash="dot"),showlegend=false),
)

pvd2

##

savefig(pvd2,"pV_and_strength.pdf",width=pw,height=ph)