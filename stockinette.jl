#=
stockinette simulation
zero-force simulation
st the box dimensions are variable
S.E. Gonzalez
written May 22 2023
updated May 24 2023
=#

using BSplines, Plots, LinearAlgebra, CSV, DataFrames


#=

Material Constants

=#
B = 0.1
k = 5
p = 3.2
rcore = 0.02
ryarn = 0.1
stitchlength = 12
targetlength = stitchlength/2

# this is the penalty for the length constraint
penalty = 10

#=

Import the initial configuration

=#
# make the b spline basis st the order is 5 and it has C3 continuity
# the second argument determines the range of t
# ALL THE SPLINES HAVE THE SAME BASIS
basis = BSplineBasis(5, 0:5)

# import the data to make the curve
initialcurve = CSV.read("filename.csv")

# make an initial spline to represent the imported curve
spl = approximate() # okay so I can't get this to work rn so maybe just import the initial control points?
# alternatively, could use the fitting control points fnc from basicbspline.jl package BUT IT NEEDS A FUNCTION??



#=

Math Functions

=#
function distance(pt1, pt2)
    d = ((pt1[1]-pt2[1])^2 + (pt1[2]-pt2[2])^2 + (pt1[3]-pt2[3])^2)^(1/2)
    return d
end

function norm(vec)
    norm = (vec[1]^2 + vec[2]^2 + vec[3]^2)^(1/2)
    return norm
end

#=

Stitch Manipulation Functions

=#

# translation functions
# INTPUT LIST IS THE X,Y,Z,T VECTOR
function translatex(list, width, cell)
    newlist = [[list[i][1] + cell*width, list[i][2], list[i][3], list[i][4]] for i in 1:length(list)]
    return newlist
end

function translatey(list, height, cell)
    newlist = [[list[i][1], list[i][2] + cell*height, list[i][3], list[i][4]] for i in 1:length(list)]
    return newlist
end

# change a knit to a purl
function mirror(list)
    newlist = [[list[i][1], list[i][2], -1*list[i][3], list[i][4]] for i in 1:length(list)]
    return newlist
end

#makes a full stitch based on just the left or right half
function fullstitch(list)
    newlist = [[-1*list[i][1], list[i][2], list[i][3], list[i][4]] for i in 1:length(list)]
    newlist = reverse(newlist) # invertes the new half so the final list moves sequentially along the stitch
    bothlists = vcat(list, newlist) #concatentates the right and left halves of the stitch
    return bothlists
end

# makes the right half of the stitch based on the left-half
function righthalf(list)
    newlist = [[-1*list[i][1], list[i][2], list[i][3], list[i][4]] for i in 1:length(list)]
    newlist = reverse(newlist) # invertes the new half so the final list moves sequentially along the stitch
    return newlist
end


#=

Bending

=#

# function to calculate bending energy
function bending(dcurve, ddcurve, dbasis, ddbasis, deltat)
    # bendingenergy = B*sum(k^2)/2 over curve

    bending = []
    fordEbend = []

    for i in 1:length(dcurve)
        # derivative of the curve crossed with the second derivative
        # returns a vector
        dcrossdd = cross(dcurve[i],ddcurve[i])

        # norm of that cross product
        normcross = norm(dcrossdd)

        # norm of the derivative of the curve
        dnorm = norm(dcurve[i])

        # incremental bending energy
        bendhere = normcross^2 * dnorm^(-5)

        # multiplies by B/2 and deltat for the arc length conversion
        push!(bending, B*bendhere*deltat/2)

        t = (i-1)*deltat

        # list of stuff for dE_bend/dP_i and the length constraint
        push!(fordEbend, [dcurve[i][1],dcurve[i][2],dcurve[i][3],ddcurve[i][1],ddcurve[i][2],ddcurve[i][3],dbasis[i],
        ddbasis[i],dcrossdd[1], dcrossdd[2], dcrossdd[3], dnorm, normcross, t])
    end

    bendingenergy = sum(bending)

    return bendingenergy, fordEbend
end

# change in bending energy wrt a control point
# dimension is x=1,y=2,z=3
# knot indicates the corrleted basis function
# tmin and tmax determine the scope of the control point wrt the parameterization
function dEbenddPi(lista, deltat, tmin, tmax, dimension, knot)
    dEbendlist = []

    list = []


    # isolate the data affected by the control point
    for i in 1:length(lista)
        if tmin <= lista[i][14] <= tmax
            push!(list, lista[i])
        end
    end

    # assumes all splines have the same basis functions

    #firstterm is d/dP_i of the norm of the cross product of the derivative and second derivative of the curve
    #second term is d/dP_i of the norm of the derivative^(-5)
    # all in accordance with the derivations I typed up and sent to you (Eq 33-35)

    #= 
    lista = [dcurve[i][1],dcurve[i][2],dcurve[i][3],ddcurve[i][1],ddcurve[i][2],ddcurve[i][3],dbasis[i],
    ddbasis[i],dcrossdd[1], dcrossdd[2], dcrossdd[3], dnorm, normcross, t] 
    =#

    for i in 1:length(list)
        if dimension == 1
            firstterm = list[i][7][knot]*(list[i][11]*list[i][5]-list[i][10]*list[i][6]) + list[i][8][knot]*(list[i][10]*list[i][3] - list[i][11]*list[i][2])
            secondtermnum = list[i][1]*list[i][7][knot]
        elseif dimension == 2
            firstterm = list[i][7][knot]*(list[i][9]*list[i][6]-list[i][11]*list[i][4]) + list[i][8][knot]*(list[i][11]*list[i][1] - list[i][9]*list[i][3])
            secondtermnum = list[i][2]*list[i][7][knot]
        elseif dimension == 3
            firstterm = list[i][7][knot]*(list[i][10]*list[i][4]-list[i][9]*list[i][5]) + list[i][8][knot]*(list[i][9]*list[i][2] - list[i][10]*list[i][1])
            secondtermnum = list[i][3]*list[i][7][knot]
        end

        secondterm = list[i][13]^2*secondtermnum/(list[i][12]^2)

        dEbend = B*deltat*(2*firstterm - 5*secondterm)/(2*list[i][12]^5)
        push!(dEbendlist, dEbend)
    end

    dEbenddP = sum(dEbendlist)
    return dEbenddP
end

#=

Compression

=#
# function to calculate force and energy densities from compression
function floof(dist)
    outer_dia = 2*ryarn
    inner_dia = 2*rcore
    dia_diff = outer_dia - inner_dia
    zeta = (dist - inner_dia)/dia_diff
    EPS = 1*10^(-8)

    pot(x) = k * dia_diff^2 * (x^(1-p) - 1 - (p - 1) * (1 - x)) / (p*(p-1))
    force(x) = k * dia_diff * (x^(-p) - 1) / p

    if zeta >= 1
        v = 0
        f = 0
    elseif EPS < zeta < 1
        v = pot(zeta)
        f = force(zeta)
    else
        v = pot(EPS) + 100 * ((dist - (inner_dia + EPS*dia_diff))/(inner_dia + EPS*dia_diff))^2
        f = force(zeta) + 100 * ((dist - (inner_dia + EPS*dia_diff))/(inner_dia + EPS*dia_diff))^2
    end

    return [v,f]
end


# THIS ONE NEEDS MASSIVE CHANGES
# so with stockinette, where will the contacts be?
# let's pretend you are only sweeping over one half of the stitch
# you have self-contacts? i don't think so
# you have right-left contacts YES
# you have up/down 1 contacts YES
# you have up/down 2 contacts YES
# you have left/right 1 contacts YES

# LET'S TRY AND ONLY DO UP AND TO THE LEFT AND THEN MULTIPLY BY TWO???? 
# SO I WOULD JUST NEED TO RECORD TWO CONTACTS FOR EACH CONTACT AND MAKE SURE THE SIGNS ARE RIGHT FOR ALL THE THINGS GOING INTO DECOMP

# THIS IS THE ONLY PLACE YOU NEED TO DO THE TRANSLATE AND MIRROR STUFF
# FOR DECOMPDPI THAT WILL BE TAKEN INTO CONSIDERATION HERE
# calculates the total compression and sets up for things needed in the gradient descent
function compression(curve, dcurve, deltat, thebasis, dbasis, height, width)
    # deltat gives the number of subdivisions; 0.01 gives 100 subdivisions for each bezier curve
    totcompeng = 0

    fordEcomp = []

    # need to make list of all other points that it could be in contact with
    #right half of the stitch
    right00 = righthalf(curve)

    #the other lists fof the right halves
    rightdcurve = reverse(dcurve)
    rightthebasis = reverse(thebasis)
    rightdbasis = reverse(dbasis)

    # the stitches above
    # left indicates that it's the left side of the stitch
    left01 = translatey(curve, height, 1)
    left02 = translatey(curve, height, 2)

    # the stitch to the left
    # right indicates that it's the right side of the stitch
    right_10 = translatex(right00, width, -1)

    # find all the contacts
    for i in 1:length(curve)
        normdr1 = norm(dcurve1[i])

        # position
        r1x = curve1[i][1]
        r1y = curve1[i][2]
        r1z = curve1[i][3]

        # derivative of the curve along all points
        dr1x = dcurve1[i][1]
        dr1y = dcurve1[i][2]
        dr1z = dcurve1[i][3]

        # contacts with right side
        for j in 1:length(right00)
            normdr2 = norm(rightdcurve[j])

            # positions
            r2x = right00[j][1]
            r2y = right00[j][2]
            r2z = right00[j][3]

            # derivatives
            dr2x = rightdcurve[j][1]
            dr2y = rightdcurve[j][2]
            dr2z = rightdcurve[j][3]

            # gives the t location for each point on each curve
            t1 = curve[i][4]
            t2 = right00[j][4]

            # distance between the points on each curve
            dist = distance(curve1[i], right00[j])

            if dist <= 2.1*ryarn
                contact = floof(dist)

                # to get potential energy, multiply the energy density by an area and then the coord transform
                potentialdensity = contact[1]
                potential = potentialdensity*normdr1*normdr2*deltat^2
                
                # this forcemag is still a density
                forcemag = contact[2]

                # direction of the force
                rhat = [(r1x-r2x), (r1y-r2y), (r1z-r2z)]
                rhat = rhat./dist

                
                totcompeng += potential


                # for use in computing dE/dP_I
                # force on left side
                push!(fordEcomp, [t1,normdr1,dr1x,dr1y,dr1z,t2,normdr2,dr2x,dr2y,dr2z,forcemag,rhat,potentialdensity, 
                thebasis[i], dbasis[i], thebasis[j], dbasis[j],r1x, r1y, r1z, r2x, r2y, r2z])
                # force on right side
                push!(fordEcomp, [t1,normdr1,dr1x,dr1y,dr1z,t2,normdr2,dr2x,dr2y,dr2z,forcemag,-1*rhat,potentialdensity, 
                thebasis[i], dbasis[i], thebasis[j], dbasis[j],r1x, r1y, r1z, r2x, r2y, r2z])
            end
        end
        
        # contacts with up/down 1
        for j in 1:length(left01)
            normdr2 = norm(dcurve[j])

            # positions
            r2x = left01[j][1]
            r2y = left01[j][2]
            r2z = left01[j][3]

            # derivatives
            dr2x = dcurve[j][1]
            dr2y = dcurve[j][2]
            dr2z = dcurve[j][3]

            # gives the t location for each point on each curve
            t1 = curve[i][4]
            t2 = left01[j][4]

            # distance between the points on each curve
            dist = distance(curve1[i], left01[j])

            if dist <= 2.1*ryarn
                contact = floof(dist)

                # to get potential energy, multiply the energy density by an area and then the coord transform
                potentialdensity = contact[1]
                potential = potentialdensity*normdr1*normdr2*deltat^2
                
                # this forcemag is still a density
                forcemag = contact[2]

                # direction of the force for up
                rhat = [(r1x-r2x), (r1y-r2y), (r1z-r2z)]
                rhat = rhat./dist

                # direction of the force for down
                minusrhat = [(r1x-r2x), -1*(r1y-r2y), (r1z-r2z)]
                minusrhat = minusrhat./dist

                # adds potential for up and down
                totcompeng += 2*potential


                # for use in computing dE/dP_I
                # the 1 up one
                push!(fordEcomp, [t1,normdr1,dr1x,dr1y,dr1z,t2,normdr2,dr2x,dr2y,dr2z,forcemag,rhat,potentialdensity, 
                thebasis[i], dbasis[i], thebasis[j], dbasis[j],r1x, r1y, r1z, r2x, r2y, r2z])
                # the 1 down one
                push!(fordEcomp, [t1,normdr1,dr1x,dr1y,dr1z,t2,normdr2,dr2x,dr2y,dr2z,forcemag,minusrhat,potentialdensity, 
                thebasis[i], dbasis[i], thebasis[j], dbasis[j],r1x, r1y, r1z, r2x, r2y-2*height, r2z])
            end
        end

        # contacts with up/down 2
        for j in 1:length(left02)
            normdr2 = norm(dcurve[j])

            # positions
            r2x = left02[j][1]
            r2y = left02[j][2]
            r2z = left02[j][3]

            # derivatives
            dr2x = dcurve[j][1]
            dr2y = dcurve[j][2]
            dr2z = dcurve[j][3]

            # gives the t location for each point on each curve
            t1 = curve[i][4]
            t2 = left02[j][4]

            # distance between the points on each curve
            dist = distance(curve1[i], left02[j])

            if dist <= 2.1*ryarn
                contact = floof(dist)

                # to get potential energy, multiply the energy density by an area and then the coord transform
                potentialdensity = contact[1]
                potential = potentialdensity*normdr1*normdr2*deltat^2
                
                # this forcemag is still a density
                forcemag = contact[2]

                # direction of the force for up
                rhat = [(r1x-r2x), (r1y-r2y), (r1z-r2z)]
                rhat = rhat./dist

                # direction of the force for down
                minusrhat = [(r1x-r2x), -1*(r1y-r2y), (r1z-r2z)]
                minusrhat = minusrhat./dist

                # adds potential for up and down
                totcompeng += 2*potential

                # for use in computing dE/dP_I
                # the 2 up one
                push!(fordEcomp, [t1,normdr1,dr1x,dr1y,dr1z,t2,normdr2,dr2x,dr2y,dr2z,forcemag,rhat,potentialdensity, 
                thebasis[i], dbasis[i], thebasis[j], dbasis[j],r1x, r1y, r1z, r2x, r2y, r2z])
                # the 2 down one
                push!(fordEcomp, [t1,normdr1,dr1x,dr1y,dr1z,t2,normdr2,dr2x,dr2y,dr2z,forcemag,minusrhat,potentialdensity, 
                thebasis[i], dbasis[i], thebasis[j], dbasis[j],r1x, r1y, r1z, r2x, r2y-4*height, r2z])
            end
        end

        # contacts with left 1, right side of stitch
        for j in 1:length(right_10)
            normdr2 = norm(rightdcurve[j])

            # positions
            r2x = right_10[j][1]
            r2y = right_10[j][2]
            r2z = right_10[j][3]

            # derivatives
            dr2x = rightdcurve[j][1]
            dr2y = rightdcurve[j][2]
            dr2z = rightdcurve[j][3]

            # gives the t location for each point on each curve
            t1 = curve[i][4]
            t2 = right_10[j][4]

            # distance between the points on each curve
            dist = distance(curve1[i], right_10[j])

            if dist <= 2.1*ryarn
                contact = floof(dist)

                # to get potential energy, multiply the energy density by an area and then the coord transform
                potentialdensity = contact[1]
                potential = potentialdensity*normdr1*normdr2*deltat^2
                
                # this forcemag is still a density
                forcemag = contact[2]

                # direction of the force
                rhat = [(r1x-r2x), (r1y-r2y), (r1z-r2z)]
                rhat = rhat./dist

                totcompeng += potential

                # for use in computing dE/dP_I
                push!(fordEcomp, [t1,normdr1,dr1x,dr1y,dr1z,t2,normdr2,dr2x,dr2y,dr2z,forcemag,rhat,potentialdensity, 
                thebasis[i], dbasis[i], thebasis[j], dbasis[j],r1x, r1y, r1z, r2x, r2y, r2z])
            end
        end
    end

    return totcompeng, fordEcomp
end

# for computing the change in compression energy wrt a control point
# tmin and tmax are the bounds of the parameter affected by the control point
# dimension is x=1,y=2,z=3
# knot indicates which basis function is affected by the control point
# curvenumber is literally an identifying for which curve is being affected
# remember that the og compression fnc sweeps over curve1 then curve2


# DOES THIS CHANGE????? YES CAUSE OF THE CURVENUMBER
function dEcompdPi(list, deltat, tmin, tmax, dimension, knot)
    if length(list) == 0
        return 0
    else
        relevantdata = []

        dEcomp = []

        # picks out all the data affected by changing the control point
        for i in 1:length(list)
            if tmin <= list[i][1] < tmax && list[i][11] != 0
                push!(relevantdata, list[i])
            elseif tmin <= list[i][6] < tmax && list[i][11] != 0
                push!(relevantdata, list[i])
            end
        end

        #=
        relevantdata = [t1,normdr1,dr1x,dr1y,dr1z,t2,normdr2,dr2x,dr2y,dr2z,forcemag,rhat,potentialdensity, 
        thebasis[i], dbasis[i], thebasis[j], dbasis[j]]
        =#

        for i in 1:length(relevantdata)
            #vec1 is the difference in basis functions
            if dimension == 1
                vec1 = [relevantdata[i][14][knot]-relevantdata[i][16][knot],0,0]
                dr1 = relevantdata[i][3]
                dr2 = relevantdata[i][8]
            elseif dimension == 2
                vec1 = [0,relevantdata[i][14][knot]-relevantdata[i][16][knot],0]
                dr1 = relevantdata[i][4]
                dr2 = relevantdata[i][9]
            elseif dimension == 3
                vec1 = [0,0,relevantdata[i][14][knot]-relevantdata[i][16][knot]]
                dr1 = relevantdata[i][5]
                dr2 = relevantdata[i][10]
            end


            # Eq 51 in the derivations write-up

            #=
            relevantdata = [t1,normdr1,dr1x,dr1y,dr1z,t2,normdr2,dr2x,dr2y,dr2z,forcemag,rhat,potentialdensity, 
            thebasis[i], dbasis[i], thebasis[j], dbasis[j]]
            =#

            # d/dP_i of the potential density
            # compression force dotted with difference in basis 
            # force * normdr1 * normdr2 * (difference in basis functions dotted w/ Rhat)
            firstterm = -1*relevantdata[i][11]*relevantdata[i][2]*relevantdata[i][7]*dot(vec1,relevantdata[i][12])

            # these are the d/dP_i of the normdr1 and normdr2
            secondterm1 = dr1*relevantdata[i][15][knot]*relevantdata[i][7]/relevantdata[i][2]
            secondterm2 = dr2*relevantdata[i][17][knot]*relevantdata[i][2]/relevantdata[i][7]

            # multiply those d/dP_i norm terms with the potential density for the product rule
            secondterm = relevantdata[i][13]*(secondterm1 + secondterm2)

            # deltat is squared because you are summing over both curves
            dEdP = deltat^2*(firstterm + secondterm)

            push!(dEcomp, dEdP)
        end

        if length(dEcomp) == 0
            return 0
        else
            dEcompdPi = sum(dEcomp)
            return dEcompdPi
        end

    end

end

#=

Length

=#
# gives length of left half
function totallength(list, deltat, tmax)
    # list should be fordEbend
    # tmax is just the last value of the parameterization
    dL = []
    for i in 1:length(list)
        if list[i][14] < tmax
            # norm of the derivative of the curve
            push!(dL, list[i][12]*deltat)
        end
    end

    totlength = sum(dL)
    return totlength
end

# change in length wrt the control point
# list should be fordEbend, knot is the part of the curve changed by the control point
# dimension x=1 y=2 z=3
# tmin and tmax give the part of the parameterization affected by the control point
function dLdPi(list, deltat, tmin, tmax, dimension, knot, targetlength, currentlength)

    relevant = []

    # picks out the data affected by the control point
    # only grabs the derivatives of the curve [1-3], the derivatives of the basis [7], and the norm of the derivative of the curve [12]
    for i in 1:length(list)
        if tmin <= list[i][14] < tmax
            push!(relevant, [[list[i][1],list[i][2],list[i][3]], list[i][7], list[i][12]])
        end
    end

    dLdPlist = []

    
    for i in 1:length(relevant)
        # dcurve * dbasis * deltat / dnorm
        length = relevant[i][1][dimension]*relevant[i][2][knot]*deltat/relevant[i][3]
        push!(dLdPlist, length)
    end

    dLdP = sum(dLdPlist)

    constant = 2*(currentlength - targetlength)
    #println(constant)
    dLdPvalue = constant*dLdP
    return dLdP, dLdPvalue
end




#=

Generate the basis functions and their derivatives

=#

function constructbasis5(basis, deltat, tmin, tmax, d)
    # constructs list of basis values for all basis functions at all values of t
    # order 5
    # only has t 0:5 for nine control points
    # it essentially just pads the output of the basis function so the indexing of the array is consistent
    units = 1/deltat
    basislist = []
    if tmin <= 0 && tmax >= 1
        for t in 0:(units-1)
            autobasis = bsplines(basis, deltat*t, Derivative(d))
            push!(basislist,[autobasis[1], autobasis[2], autobasis[3], autobasis[4], autobasis[5], 0, 0, 0, 0])
        end
    end
    if tmin<=1 && tmax>=2
        for t in units:(2*units-1)
            autobasis = bsplines(basis, deltat*t, Derivative(d))
            push!(basislist,[0, autobasis[1], autobasis[2], autobasis[3], autobasis[4], autobasis[5], 0, 0, 0])
        end
    end
    if tmin<=2 && tmax>=3
        for t in 2*units:(3*units-1)
            autobasis = bsplines(basis, deltat*t, Derivative(d))
            push!(basislist,[0, 0, autobasis[1], autobasis[2], autobasis[3], autobasis[4], autobasis[5], 0, 0])
        end
    end
    if tmin<=3 && tmax>=4
        for t in 3*units:(4*units-1)
            autobasis = bsplines(basis, deltat*t, Derivative(d))
            push!(basislist,[0, 0, 0, autobasis[1], autobasis[2], autobasis[3], autobasis[4], autobasis[5], 0])
        end
    end
    if tmin<=4 && tmax>=5
        for t in 4*units:5*units
            autobasis = bsplines(basis, deltat*t, Derivative(d))
            push!(basislist,[0, 0, 0, 0, autobasis[1], autobasis[2], autobasis[3], autobasis[4], autobasis[5]])
        end
    end


    return basislist
end


#=

Stitch Energy

=#

# gives energy of the left half of the stitch
function totalenergy(cpt1, cpt2, cpt3)
    # contructs the splines
    spline1 = Spline(basis, cpt1)
    spline2 = Spline(basis, cpt2)
    spline3 = Spline(basis, cpt3)

    # gets values for the curve positions, and 1st 2nd derivatives
    curve = [[spline1(t), spline2(t), spline3(t)] for t in 0:deltat:5]
    dcurve = [[spline1(t, Derivative(1)), spline2(t, Derivative(1)), 0] for t in 0:deltat:5]
    ddcurve = [[spline1(t, Derivative(2)), spline2(t, Derivative(2)), 0] for t in 0:deltat:5]

    # gets values of the basis functions and their derivatives at every point
    thebasis = constructbasis5(basis, deltat, 0, 5, 0)
    dbasis = constructbasis5(basis, deltat, 0, 5, 1)
    ddbasis = constructbasis5(basis, deltat, 0, 5, 2)

    #computes bending energy for curve
    bending = bending(dcurve, ddcurve, dbasis, ddbasis, deltat)
    # println("bending")
    # println(bending[1])

    # computes compression energy
    contact = compression(curve, dcurve, deltat, thebasis, dbasis, height, width)
    #println("contact")
    #println(contact[1])

    energy = contact[1] + bending[1]

    return energy, bending[2], contact[2]
end

# just prints the components of the energy for after the sim has run
# GIVES ENERGY OF TOTAL STITCH, NOT JUST HALF OF IT
function printtotalenergy(cpt1, cpt2, cpt3)
    # contructs the splines
    spline1 = Spline(basis, cpt1)
    spline2 = Spline(basis, cpt2)
    spline3 = Spline(basis, cpt3)

    # gets values for the curve positions, and 1st 2nd derivatives
    curve = [[spline1(t), spline2(t), spline3(t)] for t in 0:deltat:5]
    dcurve = [[spline1(t, Derivative(1)), spline2(t, Derivative(1)), 0] for t in 0:deltat:5]
    ddcurve = [[spline1(t, Derivative(2)), spline2(t, Derivative(2)), 0] for t in 0:deltat:5]

    # gets values of the basis functions and their derivatives at every point
    # need to fix
    thebasis = constructbasis5(basis, deltat, 0, 5, 0)
    dbasis = constructdbasis5(basis, deltat, 0, 5, 1)
    ddbasis = constructddbasis5(basis, deltat, 0, 5, 2)

    #computes bending energy for curve
    # muplitply by two cause it's only doing the left half
    bending = 2*bending(dcurve, ddcurve, dbasis, ddbasis, deltat)
    println("bending")
    println(bending[1])

    # computes compression energy
    contact = 2*compression(curve, dcurve, deltat, thebasis, dbasis, height, width)
    println("contact")
    println(contact[1])

    energy = contact[1] + bending[1]
    println(energy)
end


#=

Descent

=#

# HERE IS WHERE THE DESCENT IS DEFINED
# need to do something for the height
function descent(listbend, listcomp, deltat, tmin, tmax, dimension, knot, lambda,targetlength, length)
    benddescent = dEbenddPi(listbend, deltat, tmin, tmax, dimension, knot)
    compdescent = dEcompdPi(listcomp, deltat, tmin, tmax, dimension, knot)
    lengthstuff = dLdPi(listbend, deltat, tmin, tmax, dimension, knot, targetlength, length)
    lengthascent = lengthstuff[1]
    lengthpenalty = lengthstuff[2]

    # height stuff??
    # need these for all the knots?
    heightdescent = dEcompdPi(listcomp, deltat, tmin, tmax, 2, knot)


    descent =  compdescent + lambda*lengthascent +  benddescent + penalty*lengthpenalty
    return descent
end


# need to update/fix
function gradientdescent(cpt1, cpt2, cpt3, cpt4, learn_rate, conv_threshold, max_iter)

    # gets all the data for the original input curves
    oldenergy = totalenergy(cpt1,cpt2,cpt3)
    fordEbend = oldenergy[2]
    fordEcomp = oldenergy[3]
    oldenergy = oldenergy[1]

    oglength = totallength(fordEbend,0.01,5)

    #=

    #[t1,normdr1,dr1x,dr1y,dr1z,t2,normdr2,dr2x,dr2y,dr2z,forcemag,rhat,potentialdensity, 
    #thebasis[i], dbasis[i], thebasis[j], dbasis[j],r1x, r1y, r1z, r2x, r2y, r2z]

    df = DataFrame(r1x=Float64[], r1y=Float64[], r1z=Float64[],r2x=Float64[], r2y=Float64[], r2z=Float64[],rhatx=Float64[], rhaty=Float64[], rhatz=Float64[])

    for i in 1:length(fordEcomp)
        push!(df, (fordEcomp[i][18],fordEcomp[i][19],fordEcomp[i][20],fordEcomp[i][21],fordEcomp[i][22],fordEcomp[i][23],
        fordEcomp[i][12][1],fordEcomp[i][12][2],fordEcomp[i][12][3]))
    end

    #println(df)
    CSV.write("forcedirectioniteration0.csv",  df, header=false)

    #CSV.write("forcedirectioniteration" * string(iterations) * ".csv",  DataFrame(forexport), header=false)

    =#

    println("Here is the original energy")
    println(oldenergy)

    # prints the components of the energy
    printtotalenergy(cpt1, cpt2, cpt3)

    println("here is the og length of curve")
    println(oglength)

    # these are guesses
    lambda1 = 0.001
    lambda2 = 0.001

    converged = false
    iterations = 0

    while converged == false
        #=
        so P_i affects t_i < t < t_{i+k} and the knot i but when indexing it starts at one so actually knot = i+1
        endpoint of curve are given by the first and last control points
        for k=3
        descent(listbend, listcomp, deltat, tmin, tmax, dimension, knot, lambda, curvenumber)
        =#

        length = totallength(fordEbend,0.01,5)

        # control points for x-direction of curve
        cpt101 = cpt1[1] -learn_rate*descent(fordEbend, fordEcomp, deltat, 0, 1, 1, 1, lambda, 1, targetlength, length)
        cpt102 = cpt1[2] -learn_rate*descent(fordEbend, fordEcomp, deltat, 0, 2, 1, 2, lambda, 1, targetlength, length)
        cpt103 = cpt1[3] -learn_rate*descent(fordEbend, fordEcomp, deltat, 0, 3, 1, 3, lambda, 1, targetlength, length)
        cpt104 = cpt1[4] -learn_rate*descent(fordEbend, fordEcomp, deltat, 1, 3, 1, 4, lambda, 1, targetlength, length)
        cpt105 = cpt1[5] -learn_rate*descent(fordEbend, fordEcomp, deltat, 2, 3, 1, 5, lambda, 1, targetlength, length)
        # last control point of x should be 0 and stay 0

        cpt1 = [cpt101, cpt102, cpt103, cpt104, cpt105]

        # control points for y-direction of curve
        cpt201 = cpt2[1] -learn_rate*descent(fordEbend, fordEcomp, deltat, 0, 1, 2, 1, lambda, 1, targetlength, length)
        cpt202 = cpt2[2] -learn_rate*descent(fordEbend, fordEcomp, deltat, 0, 2, 2, 2, lambda, 1, targetlength, length)
        cpt203 = cpt2[3] -learn_rate*descent(fordEbend, fordEcomp, deltat, 0, 3, 2, 3, lambda, 1, targetlength, length)
        cpt204 = cpt2[4] -learn_rate*descent(fordEbend, fordEcomp, deltat, 1, 3, 2, 4, lambda, 1, targetlength, length)
        cpt205 = cpt2[5] -learn_rate*descent(fordEbend, fordEcomp, deltat, 2, 3, 2, 5, lambda, 1, targetlength, length)

        cpt2 = [cpt201, cpt202, cpt203, cpt204, cpt205]

        # control points for z-direction of curve
        cpt301 = cpt3[1] -learn_rate*descent(fordEbend, fordEcomp, deltat, 0, 1, 3, 1, lambda, 2, targetlength, length)
        cpt302 = cpt3[2] -learn_rate*descent(fordEbend, fordEcomp, deltat, 0, 2, 3, 2, lambda, 2, targetlength, length)
        cpt303 = cpt3[3] -learn_rate*descent(fordEbend, fordEcomp, deltat, 0, 3, 3, 3, lambda, 2, targetlength, length)
        cpt304 = cpt3[4] -learn_rate*descent(fordEbend, fordEcomp, deltat, 1, 3, 3, 4, lambda, 2, targetlength, length)
        cpt305 = cpt3[5] -learn_rate*descent(fordEbend, fordEcomp, deltat, 2, 3, 3, 5, lambda, 2, targetlength, length)
        # first control point of z should be 0 and stay 0

        cpt3 = [cpt301, cpt302, cpt303, cpt304, cpt305]

        #computes the new energies and the new lists used to determine dE/dP_i
        newenergy = totalenergy(cpt1, cpt2, cpt3)
        fordEbend = newenergy[2]
        fordEcomp = newenergy[3]
        newenergy = newenergy[1]

        println("here is the updated energy")
        println(newenergy)
        println(iterations)

        #gradient ascent part for the length constraint
        lambda = lambda + (totallength(fordEbend,0.01,5) - targetlength)

        if abs(oldenergy - newenergy) <= conv_threshold
            converged = true
            println("We converged!")
            println("here is the new energy")
            println(newenergy)
            println("I did this many iterations:")
            println(iterations)
        else
            oldenergy = newenergy
        end

        iterations +=1

        if iterations > max_iter
            converged = true
            println("We reached the max number of iterations!")
            println(newenergy)
        end

    end

    println("here is the end length of curve")
    println(totallength(fordEbend,0.01,5))

    println("Here are the new control points")

    println(cpt1)
    println(cpt2)
    println(cpt3)

    println("Here is the energy breakdown:")

    printtotalenergy(cpt1, cpt2, cpt3)

    return cpt1, cpt2, cpt3
end


#=

Where the code actually runs

=#

# actually doing stuff now

deltat = 0.01

# the control point lists were defined at the very beginning
# the last number is the max number of iterations
grad = gradientdescent(xpoints1, ypoints1, ypoints2, zpoints2, 0.001, 0.000000001, 100000)



#=

Visualizations

=#

#Plot the finalized clasp
xspline = Spline(basis, grad[1])
yspline = Spline(basis, grad[2])
zspline = Spline(basis, grad[3])
curve = [[xspline(t) yspline(t) zspline(t)] for t in 0:0.01:5]
plotcurve = plot(Tuple.(curve), xlabel="X", ylabel="Y", zlabel="Z")
png(plotcurve, "finalcurve.png")