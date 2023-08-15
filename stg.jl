
## A simple simulator of the syntetic turbulence generator 
## by Shur et al., Flow Turbulence and Combustion (2004) 

## There are few things, which are not implemented: 
## 1) function to find the wave number corresponding to the spectral maximum
## 2) function to find the max local scale of the most energy containg eddies
## 3) not sure that mapping of the wave space on the spatialtemporal Fourier space is fully correct - TBD 


using PyPlot

function ellipsoid(center, rx, ry, rz, ngrid=25)

    # Radii corresponding to the coefficients:
    #rx, ry, rz = 1 ./ sqrt.(coefs)
    # rx = 0.5;
    # ry = 0.6;
    # rz = 4.2;

    # Set of all spherical angles:
    u = range(0, 2pi, length=ngrid)
    v = range(0, pi, length=ngrid)

    # Cartesian coordinates that correspond to the spherical angles:
    # (this is the equation of an ellipsoid):
    x = [rx * x * y for (x, y) in  Iterators.product(cos.(u), sin.(v))]
    y = [ry * x * y for (x, y) in Iterators.product(sin.(u), sin.(v))]
    z = [rz * x * y for (x, y) in Iterators.product(ones(length(u)), cos.(v))]
    return x .+ center[1], y .+ center[2], z .+ center[3]

end

function setupWaveNumberSpace(alpha::Float64, k_min::Float64, N::Int64)::Array{Float64,1} 

    ## N -  number of modes 
    ## alpha - the constant, 0.01 <= alpha <= 0.05
    ## k_min - the minimum wavenumber in the set 
    ## k_min = 0.5*ke_min, 
    ## ke_min = 2pi/le_max,
    ## le_max - the max local scale of the most energy containg eddies over the inflow boundary

    k = [] ## wavenumber vector
    
    k0 = k_min*(1+alpha)^(1-0) 
    push!(k,k0)

    for n = 1:N
        kn = k_min*(1+alpha)^(n-1)
        push!(k,kn)
    end

    return k;
end

function f_eta(k::Array{Float64,1},k_eta::Float64)::Array{Float64,1} 
     ## damping function of the energy spectrum in the vicinity of the wave number 
     ## corresponding to the Kolmogorov length scale (\eta)
     ## k_eta = 2pi/eta
     return exp.( -1.0.*(12.0 .* k ./k_eta).^2 )
end

function f_cut(k::Array{Float64,1}, lcut::Float64)::Array{Float64,1}
    ## damping function of the energy spectrum at waves numbers 
    ## larger than the Nyquist value, k_cut
    ## assumptions: 
    ## lcut is the min grid length scale 
    kcut = 2.0*pi/lcut
    return exp.(-( 4.0 .* max.(k .- 0.9.* kcut,0.0)./kcut ).^3.0)
end


function E(k::Array{Float64,2}, ke::Float64, k_eta::Float64, l_cut::Float64)::Array{Float64,1}

    ## compute a modified von Karman spectrum 
    ## ke coresponds to the wavelength of the most energy containg mode 
    ## lcut is the min grid length scale 
    ## k_eta - wavenumber corresponding to the Kolmogorov length scale

    return ( (k./ke).^4 ./ (1.0 .+ 2.4.*(k./ke).^2).^(17.0/6.0) ) .*f_eta(k,k_eta) .*f_cut(k,l_cut)
end

function calcAij(Rij)::Array{Float64,2}

    ## a Cholesky decomposition of the Reynolds stress tensor, Rij
    
    Aij = zeros(Float64,3,3)
    Aij[1,1] = sqrt(Rij[1,1])
    Aij[2,1] = Rij[2,1]/Aij[1,1]
    Aij[2,2] = sqrt(Rij[2,2]*Rij[2,2] - Aij[2,1]*Aij[2,1])
    Aij[3,1] = Rij[3,1]/Aij[1,1] 
    Aij[3,2] = (Rij[3,2]-Aij[2,1]* Aij[3,1])/Aij[2,2]
    Aij[3,3] = sqrt(Rij[3,3] -Aij[3,1]*Aij[3,1] - Aij[3,2]*Aij[3,2]  )

    return Aij;

end

function prime()


    ## flow parameters 
    TInf = 300.0 ## temperature [K]
    PInf = 101325.0 ## pressure [Pa]
    R = 287.0 ## gas const
    rhoInf = PInf/R/TInf ## density [kg/m3]
    muInf = 1.79e-5 ## dynamic viscosity [Pa*s]
    nuInf = muInf/rhoInf ## kinematic viscosity [m2/s]
    DInf = 0.1 ## length scale [m]
    UInf = 4.0  ## bulk velcoity [m/s]
    TiInf = 0.05 ## Turbulence intensity, u'/U [-]
    uFluctInf = TiInf*UInf ## fluctuation velocity [m]
    display("Ufluct: "*string(uFluctInf)*" [m/s]")
    lFluctInf = DInf ## the length scale [m]
    epsInf = uFluctInf^3/lFluctInf ## dissipation 
    display("Dissipation: "*string(epsInf)*" [m2/s3]" )
    etaInf = (nuInf^3/epsInf)^(1/4) ## Kolmogorov length scale [m]
    display("Kolmogorov length scale: "*string(etaInf))

    ## Reynolds stresses 
    RijInf = [
        uFluctInf*uFluctInf 0.0 0.0 
        0.0  uFluctInf*uFluctInf 0.0
        0.0  0.0 uFluctInf*uFluctInf                 
    ]
    
    ## a Cholesky decomposition of the Reynolds stress
    AijInf = calcAij(RijInf)

    le_max = 25 ## TBD 
    ke_min = 2.0*pi/le_max
    k_min = 0.5*ke_min
    l_cut = 0.001 ## min grid length scale
    l_eta = etaInf  ## kolmogorov length scale
    k_eta = 2.0*pi/l_eta 
    ke = 1.0 # is the wave number corresponding to the spectral maximum , TBD 

    ## vecocity vector 
    Ux = UInf
    Uy = 0.0
    Uz = 0.0

    t = [] ## time vector
    push!(t,0)

    dt = 0.5 ## time step 

    ## arbitral point on the inflow boundary 
    Point_x = 0.0
    Point_y = 0.0
    Point_z = 0.0

    N = 1000 ## number of nodes 

    k = setupWaveNumberSpace(0.01,  k_min, N )

    Ek = E(k, ke, k_eta, l_cut)

    deltas = diff(k)

    EkInt = sum( Ek[i] * deltas[i] for i in 1:length(deltas) ) 
    #display(EkInt)

    ## compute normalized amplitudes of the modes 
    qn = Ek[1:end-1] .* deltas[1:end] ./EkInt    
    nN = size(qn,1)



    ex0, ey0, ez0 = ellipsoid([0 0 0],1.0,1.0,1.0) ## a unit sphere 
    ## defines by the radial distance, polar angle, and azimuthal angle
    ## radial distance: r ≥ 0,
    ## polar angle: 0° ≤ θ ≤ 180° (π rad) -> theta^n
    ## azimuth : 0° ≤ φ < 360° (2π rad) -> phi^n

    ## compute rand points on the surfphase of a unit sphere 
    theta = rand(Float64, nN ).*pi
    phi =   rand(Float64, nN ) .*2.0 .*pi
    radiu = ones(Float64, nN )  

    ## compute normals
    sigma_x  = radiu .* sin.(theta) .*cos.(phi)
    sigma_y  = radiu .* sin.(theta) .*sin.(phi)
    sigma_z  =  radiu .* cos.(theta) 

    kk = k[1:end-1]
    Ekk = Ek[1:end-1]
    
    ## pseudo-position vector
    rxStar = 2.0 .* pi ./ kk ./ le_max .* (Point_x .- Ux.*t)
    ryStar = Point_y 
    rzStar = Point_z 


    ## compute a syntetic spectrum 
    EkFx = Ekk./EkInt .* sigma_x .*cos.( sigma_x .* kk .* rxStar .+ phi )
    
    Uf = []
    push!(Uf,UxFluct )

    ## compute 200 time steps of the solution 

    for i = 1:200
        tau = t[end] + dt;
        beta = 2.0*sqrt(3.0/2.0)
        rxStar = 2.0 .* pi ./ kk ./ le_max .* (Point_x .- Ux.*tau)
        #ryStar = Point_y 
        #rzStar = Point_z 

        beta = 2.0*sqrt(3.0/2.0)
        Uxrms = beta*sum(sqrt.(abs.(qn)).*sigma_x .*cos.( sigma_x .* kk .* rxStar .+ phi ))
        
       # display(Uxrms)

        push!(t,tau)
        push!(Uf, Uxrms )

    end



    lw1 = 1.0;

    figure(1)
    clf()
    mesh(ex0, ey0, ez0, alpha=0.50, linewidth = lw1, ec = "g")
    plot3D(xP,yP,zP,"or",markersize = 2.0)
    xlabel("x")
    ylabel("y")
    zlabel("z")

    figure(2)
    clf()

    subplot(1,2,1)
    
    loglog(kk,EkFx,"-r",linewidth = lw1, label = "syntetic")
    loglog(k,Ek,"-g",linewidth = lw1+1.0,label = "von Karman")
    #loglog(kk,EkFy,"--b",linewidth = lw1)
    #loglog(kk,EkFz,"-.c",linewidth = lw1)
    xlabel("k")
    ylabel("E(k)")
    legend()
    grid()

    subplot(1,2,2)
    plot(t,Uf.*AijInf[1,1],"-g")
    xlabel("t")
    ylabel("Ux'")
    grid()


    figure(3)
    clf()
    subplot(2,1,1)
    loglog(k,f_cut(k,l_cut),"-g",label="f_cut")
    legend()
    subplot(2,1,2)
    loglog(k,f_eta(k,k_eta),"-g",label ="f_eta")
    legend()


end


prime()
