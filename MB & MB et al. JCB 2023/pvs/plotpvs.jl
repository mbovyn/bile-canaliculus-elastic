using CSV, DataFrames, PlotlyJS, DataInterpolations

cd(@__DIR__)
nfiles=12

scl=3.6

## get p-v curves

global sss=String[]
for i=0:nfiles-1
    push!(sss, string("p-v/","v-",i,".txt"))
end

dfs = DataFrame.(CSV.File.(sss))

##

dfsall=dfs[1]
for df in dfs[2:end]
    global dfsall=vcat(dfsall,df)
end

## get vcrit values

vcrits=DataFrame(CSV.File.("../max_expansion_w_cut_depth/vcrit.csv"))

## non-annotated p-v curves

p=plot(Layout(width=500,height=500))

for df in dfs
    addtraces!(p,scatter(df,x=:p,y=:v))
end
p

##

nobinterp=LinearInterpolation(dfs[1].p,dfs[1].v)

p0i=nobinterp.(vcrits.vcrit)

pcrit=zeros(length(sss)-1)
for (i,df) in enumerate(dfs[2:end])
    interp=LinearInterpolation(df.p,df.v)
    pcrit[i]=interp(vcrits[vcrits.c .== df[1,:bulkheadyoungs],:vcrit][1])
end

## annotated p-v curves

p=plot(Layout(width=600,height=500,
    legend_title_text="Bulkhead<br>modulus",
    xaxis_title="Pressure",
    yaxis_title="Volume"
))

icrit=3

for (i,df) in enumerate(dfs)
    addtraces!(p,scatter(df,x=:p,y=:v,name=string(df[1,:bulkheadyoungs])))
    if i>1
        addtraces!(p,scatter(y=vcrits.vcrit[i-1].*[1,1],x=[p0i[i-1],df[icrit,:p]],line=attr(color="black"),showlegend=false))
    end
end
p

#savefig(p,"pV_curves.pdf",width=600,height=500);

##

cs=unique(dfsall,:bulkheadyoungs)[!,:bulkheadyoungs]

pinc=pcrit./p0i
push!(pinc,1)
pinc=sort(pinc)

pincinterp=LinearInterpolation(pinc,cs)

## pressure increase with modulus

cmin=1.21293
cmax=10

pp=plot(Layout(
    width=350,
    height=300,
    xaxis_title="Bulkhead modulus / Cortex modulus",
    yaxis_title="Increase in pressure held<br>(constant volume)",
    yaxis_range=[1,1.61],
    xaxis_range=[0,10.3],
    yaxis=attr(tickmode="auto",nticks=7,showgrid=false,gridcolor="Grey",showline=true,linecolor="Black",ticks="inside"),
    xaxis=attr(tickmode="auto",nticks=14,showgrid=false,gridcolor="Black",showline=true,linecolor="Black",ticks="inside"),
    plot_bgcolor="White",
    #font=attr(size=14),
))
addtraces!(pp,scatter(
    x=[cmin,cmin,0,0,cmax,cmax],
    y=[1,pincinterp(cmin),pincinterp(cmin),pincinterp(cmax),pincinterp(cmax),1],
    mode="lines",fill="toself",fillcolor="Grey",showlegend=false,opacity=.2,line_color=:transparent)
)
addtraces!(pp,scatter(x=cs,y=pinc,line_color="RoyalBlue",showlegend=false))
pp

#savefig(pp,"p_inc_w_mod.pdf",width=350,height=300);

##

pvp=plot(Layout(
    xaxis_title="Volume (μm³)",
    yaxis_title="Pressure/Cortex Modulus",
    width=400,
    height=400,
    legend=attr(x=0.025, y=.975),
    xaxis=attr(tickmode="auto",nticks=10)
))

#pvp

fac=1.22
facs=1.05:.01:1.5

vpinterp1=LinearInterpolation(dfs[3].p,dfs[3].v)
p1=vpinterp1(vcrits.vcrit[2]*fac)
p1-vpinterp1(vcrits.vcrit[2])

vpinterp0=LinearInterpolation(dfs[1].p,dfs[1].v)
p0=vpinterp0(vcrits.vcrit[2]*fac)
p0-vpinterp0(vcrits.vcrit[2])

(p1-vpinterp1(vcrits.vcrit[2]))/(p0-vpinterp0(vcrits.vcrit[2]))
pmore=(vpinterp1.(vcrits.vcrit[2].*facs).-vpinterp1(vcrits.vcrit[2]))./(vpinterp0.(vcrits.vcrit[2].*facs).-vpinterp0(vcrits.vcrit[2]))

addtraces!(pvp,
scatter(
    x=[vcrits.vcrit[2],vcrits.vcrit[2]*fac,6,6]*scl^3,
    y=[pcrit[2],p1,p1,pcrit[2]],
    fill="toself",mode="lines",fillcolor="Red",showlegend=false,opacity=.2,line_color=:transparent),
scatter(
    x=[vcrits.vcrit[2],vcrits.vcrit[2]*fac,6,6]*scl^3,
    y=[vpinterp0(vcrits.vcrit[2]),p0,p0,vpinterp0(vcrits.vcrit[2])],
    fill="toself",mode="lines",fillcolor="RoyalBlue",showlegend=false,opacity=.2,line_color=:transparent),
)

addtraces!(pvp,
    scatter(dfs[3],x=dfs[3].v*scl^3,y=:p,name="With Bulkheads",line_color="Red"),
    scatter(dfs[1],x=dfs[1].v*scl^3,y=:p,name="No Bulkheads",line_color="RoyalBlue"),
)
pvp

##

plot(scatter(x=facs,y=pmore),
    Layout(yaxis_range=[1,1.25],
    xaxis_range=[1,1.55],
    xaxis_title="Increase in volume before rupture",
    yaxis_title="Δp with bulkheads / Δp without",
    width=400,
    height=300
    )
)

## find vcrits

pvc=plot(
    #Layout(xaxis_range=[3.9,4.1],yaxis_range=[.99,1.01])
)
vcrit2=zeros(length(dfs))
for (i,df) in enumerate(dfs)
    addtraces!(pvc,scatter(df,x=:v,y=df.maxx))

    vcinterp=LinearInterpolation(df.v,df.maxx)
    vcrit2[i]=vcinterp(1.25)
    addtraces!(pvc,scatter(x=[vcrit2[i]],y=[1.25],line_color="Black",showlegend=false))
end
pvc

##

pvc2=plot(
    #Layout(xaxis_range=[3.9,4.1],yaxis_range=[.99,1.01])
)
for (i,df) in enumerate(dfs)
    addtraces!(pvc2,scatter(df,x=:v,y=df.maxz))
end
pvc2

[pvc;pvc2]

## 

pxs=plot(Layout(
    xaxis_range=[.95,1.55]
))
for df in dfs
    addtraces!(pxs,scatter(df,x=:maxx,y=:p))
end
pxs

## pressure increase to Dbleb with modulus

fac=1.21
icrit=3
xcrit0=dfs[1].maxx[icrit]

xpinterp0=LinearInterpolation(dfs[1].p,dfs[1].maxx)
p0f=xpinterp0(xcrit0*fac)
p00=dfs[1].p[icrit]
dp0=p0f-p00
addtraces!(pxs,scatter(x=[1,fac]*xcrit0,y=[p00,p0f],line_color="Black"))

pmore=zeros(length(dfs[2:end]))

for (i,df) in enumerate(dfs[2:end])

    xcrit=df.maxx[icrit]

    xpinterp1=LinearInterpolation(df.p,df.maxx)
    p1f=xpinterp1(xcrit*fac)
    p10=df.p[icrit]
    dp1=p1f-p10

    pmore[i]=dp1/dp0

    addtraces!(pxs,scatter(x=[1,fac]*xcrit,y=[p10,p1f],line_color="Black"))
end

##

sinterp=LinearInterpolation(pmore,cs[2:end])
satc=sinterp(cmin)

pw=400
ph=350

swp=plot([scatter(x=[cmin,cmin,0],y=[1,satc,satc],text=["","",string(round(satc,digits=2))],textposition="left",showlegend=false,line_color="Grey",mode="lines"),
    scatter(x=cs[2:end],y=pmore,line_color="RoyalBlue",mode="markers",showlegend=false),
    scatter(x=cs[2]:.01:cs[end],y=sinterp.(cs[2]:.01:cs[end]),line_color="RoyalBlue",showlegend=false)],
    Layout(
    xaxis_title="Bulkhead modulus / Cortex modulus",
    yaxis_title="Δp / Δp₀",
    width=pw,
    height=ph,
    #xaxis_range=[0,10.5],
    xaxis_range=[-1,1.1],
    yaxis_range=[1,2.4],
    xaxis_type="log",
    xaxis=attr(tickmode="auto",nticks=7),
    yaxis=attr(tickmode="auto",nticks=10),
    )
)

##
#savefig(swp,"s_with_c.pdf",width=pw,height=ph)

## pressure increase above p at dstar

pscl=400/300
pw=ceil(Int64,200*pscl)
ph=ceil(Int64,250*pscl)

pstar0=dfs[1][icrit,:p]

pstars=zeros(length(dfs)-1)
for (i,df) in enumerate(dfs[2:end])
    pstars[i]=df[icrit,:p]
end

dp=pstars./pstar0

dpinterp=LinearInterpolation(dp,cs[2:end])
dpp=dpinterp(cmin)

strwc=plot([scatter(x=[cmin,cmin,0],y=[1,dpp,dpp],showlegend=false,line_color="Grey",mode="lines"),
scatter(x=cs[2:end],y=dp,line_color="RoyalBlue",mode="markers",showlegend=false),
scatter(x=cs[2]:.01:cs[end],y=dpinterp.(cs[2]:.01:cs[end]),line_color="RoyalBlue",showlegend=false)],
Layout(
    xaxis_title="Bulkhead modulus / Cortex modulus",
    yaxis_title="Δp",
    width=pw,
    height=ph,
    #xaxis_range=[0,10.5],
    xaxis_range=[-1,1.1],
    yaxis_range=[1,5],
    xaxis_type="log",
    xaxis=attr(tickmode="auto",nticks=7),
    yaxis=attr(tickmode="auto",nticks=10),
    font=attr(size=8*pscl,family="Arial")
    )
)
strwc

##

savefig(strwc,"strength_with_c.pdf",width=pw,height=ph)