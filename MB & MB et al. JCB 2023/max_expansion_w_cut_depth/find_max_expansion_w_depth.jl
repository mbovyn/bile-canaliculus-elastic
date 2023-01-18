using CSV, DataFrames, PlotlyJS, DataInterpolations

cd(@__DIR__)
nfiles=10

##

pres=String[]
push!(pres,string("inflated_sweep_max"))
push!(pres,string("inflated_sweep_med_low"))
push!(pres,string("inflated_sweep_med"))
push!(pres,string("inflated_sweep_med_high"))
push!(pres,string("inflated_sweep_min"))
scl=3.6
ds=(.25 .+ (1.18-.25).*[1,3/4,1/2,1/4,0])/scl

dfs=Array[]
dfse=Array[]
df = [DataFrame() for _ in 1:5]

for s=1:5
    # get original data

    global sss=String[]
    for i=0:nfiles-1
        push!(sss, string(pres[s],"/v-",i,".txt"))
    end

    push!(dfs,DataFrame.(CSV.File.(sss)))

    # construct a single dataframe with all original data
    #push!(df,vcat(dfs[s][1],dfs[s][2]))
    df[s]=vcat(dfs[s][1],dfs[s][2])
    for i=3:length(sss)
        df[s]=vcat(df[s],dfs[s][i])
    end

    ## get expansion data

    global ssse=String[]
    for i=0:nfiles-1
        push!(ssse,string(pres[s],"/v-",i,"_expand.txt"))
    end

    push!(dfse,DataFrame.(CSV.File.(ssse)))

end

## explore data

# plot(scatter(x=df[(df.bulkheadyoungs .== .01),:][:,:v],y=df[(df.bulkheadyoungs .== .01),:][:,:maxx]))
# gbby=groupby(df,:bulkheadyoungs)
# p=plot(scatter(gbby[1],x=:v,y=:maxx))
# for i=2:length(gbby)
#     addtraces!(p,scatter(gbby[i],x=:v,y=:maxx))
# end
# p

## plot vertices of outside of bulkhead before and after ablation

# which C
i=4
#which V
j=4
#which d
k=3

vs=unique(dfse[k][i],:v)[:,:v]
ebv=groupby(sort(unique(dfse[k][i]),:x),[:v,:ablated])
plot([scatter3d(ebv[(vs[j],0)],x=:x,y=:y,z=:z,name="before"),
    scatter3d(ebv[(vs[j],1)],x=:x,y=:y,z=:z,name="after"),
    scatter3d(
        x=ebv[(vs[j],0)][1,:nickc].+[-1,0,1]*ebv[(vs[j],0)][1,:nickw]/2,
        y=ebv[(vs[j],0)][1,:y].*[1,1,1],
        z=[ebv[(vs[j],0)][1,:bbmin],ebv[(vs[j],0)][1,:nickdepth],ebv[(vs[j],0)][1,:bbmin]],
        showlegend=false),
    scatter3d(
        x=ebv[(vs[j],0)][1,:nickc].+[-1,1]*.75,
        y=ebv[(vs[j],0)][1,:y].*[1,1],
        z=ebv[(vs[j],0)][1,:bbmin].*[1,1],
        showlegend=false,mode="lines",line=attr(width=4,color="black",dash="dash"))
    ],
    Layout(
    xaxis_range=[minimum(ebv[(vs[j],0)].x),maximum(ebv[(vs[j],0)].x)],
    scene_aspectratio=attr(x=maximum(ebv[(vs[j],0)].x)-minimum(ebv[(vs[j],0)].x),
        y=maximum(ebv[(vs[j],0)].y)-minimum(ebv[(vs[j],0)].y),
        z=maximum([maximum(ebv[(vs[j],0)].z),ebv[(vs[j],0)][1,:nickdepth]])-minimum(ebv[(vs[j],0)].z)),
        title=string("Bulkhead young's = ",dfse[k][i][1,:bulkheadyoungs],", v = ",vs[j])
    #height=400,
    #width=500
))

##

bap=plot([scatter(ebv[(vs[j],0)],x=ebv[(vs[j],0)][!,:x]*scl,y=ebv[(vs[j],0)][!,:z]*scl,name="before"),
    scatter(ebv[(vs[j],1)],x=ebv[(vs[j],1)][!,:x]*scl,y=ebv[(vs[j],1)][!,:z]*scl,name="after"),
    scatter(
        x=(ebv[(vs[j],0)][1,:nickc].+[-1,0,1]*ebv[(vs[j],0)][1,:nickw]/2)*scl,
        y=[ebv[(vs[j],0)][1,:bbmin],ebv[(vs[j],0)][1,:nickdepth],ebv[(vs[j],0)][1,:bbmin]]*scl,
        name="ablated region",fill="toself",mode="lines",line=attr(color="grey")),
    scatter(
        x=(ebv[(vs[j],0)][1,:nickc].+[-1,1]*.75)*scl,
        y=ebv[(vs[j],0)][1,:bbmin].*[1,1]*scl,
        showlegend=false,mode="lines",line=attr(color="black",dash="dash"))
    ],
    Layout(
        xaxis_range=[minimum(ebv[(vs[j],1)].x*scl),maximum(ebv[(vs[j],1)].x*scl)]+[-.1,.1],
        yaxis=attr(scaleanchor="x", scaleratio=1),
        height=400,
        width=600,
        #legend=attr(x=.05, y=.95),
        xaxis_title="x (μm)",
        yaxis_title="z (μm)",
        title=string("Bulkhead modulus / cortex modulus = ",dfse[k][i][1,:bulkheadyoungs],", Volume = ",round(Int64,vs[j]*scl^3)," μm³"),
    )
)

#savefig(bap,"bulkhead_before_after.pdf",width=600,height=400)

## function to get expansion ratios

getlr(vec)=[vec[1:floor(Int64,length(vec)/2)],vec[ceil(Int64,length(vec)/2):length(vec)]]
lrxz(x,z)=[getlr(x),getlr(z)]

function expansion_ratio(fractionin,before,after)
    lrb=lrxz(before[:,:x],before[:,:z])
    lra=lrxz(after[:,:x],after[:,:z])

    interpbl=LinearInterpolation(lrb[1][1],lrb[2][1])
    interpbr=LinearInterpolation(reverse(lrb[1][2]),reverse(lrb[2][2]))
    interpal=LinearInterpolation(lra[1][1],lra[2][1])
    interpar=LinearInterpolation(reverse(lra[1][2]),reverse(lra[2][2]))

    expansion(z)=(interpar(z)-interpal(z))/(interpbr(z)-interpbl(z))

    zlook=minimum(before[:,:z]).+fractionin.*(maximum(before[:,:z])-minimum(before[:,:z]))
    #zlook=before[1,:bbmin]
    #print(zlook)
    return expansion.(zlook)
    #return zlook
end

## plot expansion ratio as a function of depth into the bulkhead

# perf=plot()
# for j in 1:9
#     addtraces!(perf,scatter(x=0:.01:1,y=expansion_ratio(0:.01:1,ebv[(vs[j],0)],ebv[(vs[j],1)])))
# end
# perf

## create a dataframe with the expansion values, volumes and bulkhead youngs

edf = [DataFrame() for _ in 1:5]
newdf= [DataFrame() for _ in 1:5]
gg = Array[]

for s=1:5

    #push!(edf,DataFrame[])
    for (i,subdfse) in enumerate(dfse[s])

        global eb=groupby(sort(unique(subdfse),:x),[:v])

        global es=Float64[]
        global vs=Float64[]
        for sub in eb
            before=sub[sub.ablated .== 0,:]
            after=sub[sub.ablated .== 1,:]
            push!(es,expansion_ratio(.2,before,after))
            push!(vs,unique(sub,:v)[:,:v][1])
        end
        #print(es,"\n")
        global edf[s]=vcat(edf[s],DataFrame(v=vs,expansion=es,bulkheadyoungs=subdfse[1,:bulkheadyoungs]))
    end



    # create a joined dataframe including original data and expansion

    newdf[s]=innerjoin(edf[s],df[s],on = [:bulkheadyoungs,:v])

end

## plot expansion with volume

pe=plot(Layout(
    xaxis_title="Volume (internal units)",
    yaxis_title="Expansion ratio",
    legend_title_text="Bulkhead modulus"
    )
)
for subgg in gg
    addtraces!(pe,scatter(subgg,x=:v,y=:expansion,name=string(subgg.bulkheadyoungs[1])))
end
pe

## plot circularity with volume

toplot=5
gg=groupby(newdf[toplot],:bulkheadyoungs)

pc=plot(Layout(xaxis_title="Volume",
yaxis_title="width/height")
)
for subgg in gg
    addtraces!(pc,scatter(subgg,x=:v,y=subgg.maxx./subgg.maxz,name=string(subgg.bulkheadyoungs[1])))
end
pc

## find critical volumes

vcrit=Array[]
cvc=Array[]

for s=1:5

    gg=groupby(newdf[s],:bulkheadyoungs)

    push!(vcrit,Float64[])
    for subgg in gg
        interp = LinearInterpolation(subgg.v,subgg.maxx/subgg.maxz[1])
        push!(vcrit[s],interp(1))
    end
    #vcrit

    # find expansions at the critical volumes

    push!(cvc,Float64[])
    for (i,subgg) in enumerate(gg)
        interp=LinearInterpolation(subgg.expansion,subgg.v)
        push!(cvc[s],interp(vcrit[s][i]))
    end
    #cvc

end

## plot expansions with volume with critical volumes overlayed

toplot=5
gg=groupby(newdf[toplot],:bulkheadyoungs)

pe2=plot(Layout(
    xaxis_title="Volume (internal units)",
    yaxis_title="Expansion ratio",
    legend_title_text="Bulkhead modulus"
    )
)
for subgg in gg
    addtraces!(pe2,scatter(subgg,x=:v,y=:expansion,name=string(subgg.bulkheadyoungs[1])))
end
addtraces!(pe2,scatter(x=vcrit[toplot],y=cvc[toplot]))
pe2

## plot expansion at the critical volume against the bulkhead modulus

pw=500
ph=350
pf=plot(
    Layout(
        xaxis_title="Bulkhead modulus / cortex modulus",
        yaxis_title="Expansion ratio",
        yaxis_range=[1,1.4],
        #xaxis_range=[0,15.1],
        xaxis_type="log",
        #xaxis_range=[-2.1,1.1],
        #showlegend=false,
        yaxis=attr(showgrid=false,gridcolor="Grey",showline=true,linecolor="Black",ticks="inside"),
        xaxis=attr(showgrid=false,gridcolor="Black",showline=true,linecolor="Black",ticks="inside"),
        #paper_bgcolor="Blue",
        plot_bgcolor="White",
        width=pw,
        height=ph,
        legend_title_text="Cut depth (μm)"
    )
)


c=Array[]
for s=1:5
    push!(c,unique(df[s],:bulkheadyoungs)[!,:bulkheadyoungs])
    erinterp=LinearInterpolation(c[s],cvc[s])
    rmax=1.25
    rmean=1.11
    rmin=1.04
    cmin=erinterp(rmin)
    cmax=erinterp(rmax)
    cmean=erinterp(rmean)
    #addtraces!(pf,scatter(x=[0,0,cmax,cmax,cmin,cmin,0], y=[rmin,rmax,rmax,1,1,rmin,rmin], fill="toself",fillcolor="Grey",opacity=.2,line_color="Grey",mode="lines"))
    #ddtraces!(pf,scatter(x=[0,cmean,cmean],y=[rmean,rmean,1],line=attr(dash="dash",color="black")))
    addtraces!(pf,scatter(x=c[s],y=cvc[s],name=string(round(ds[s]*1.5*3.6,digits=2))))
end
#print("Bulkhead modulus corresponding to mean expansion ratio is ",cmean,"\n")
pf

##

#savefig(pf,"expansion_ratio_v_modulus.pdf",width=pw,height=ph)

##

#tosave=DataFrame(c=c,cvc=cvc)
#CSV.write("medcut.csv",tosave)

##

plot(ds,maximum.(cvc))

##
tosave=DataFrame()
for s=1:5
    global tosave=vcat(tosave,DataFrame(d=ds[s],expansion_ratio=cvc[s],modulus=c[s]))
end
#CSV.write("expansions.csv",tosave)

## save critical volumes for later

#plot(scatter(y=vcrit[1],x=c[1]),Layout(xaxis_type="log"))
pvcs=plot(scatter(y=vcrit[1][2:7],x=c[1][2:7]))
vcritinterp=LinearInterpolation(vcrit[1][2:7],c[1][2:7]);

crange=[.1,1,2,3,4,5,6,7,8,9,10]
vcritdata=DataFrame(c=crange,vcrit=vcritinterp.(crange))
addtraces!(pvcs,scatter(vcritdata,x=:c,y=:vcrit))
pvcs
##
#CSV.write("vcrit.csv",vcritdata)