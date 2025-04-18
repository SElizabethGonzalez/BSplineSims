#=
stockinette simulation
zero-force simulation
st the box dimensions are variable
S.E. Gonzalez
written May 22 2023
updated May 31 2023
=#

#=
right now, this simulation is essentially fixing the height and letting the width vary
to fix the width, implement the fixwidth fnc in the gradient descent and comment out the first x cpt
=#

using BSplines, Plots, LinearAlgebra, CSV, DataFrames


#=

Material Constants

=#
B = 0.045
k = 0.0006
p = 2.4
rcore = 0.3
ryarn = 0.74
stitchlength = 10.4
targetlength = stitchlength/2

width = 2.68
height = 1.829

# this is the penalty for the length constraint
penalty = 10

# number of subdivisions per t=1 unit for contact calculations
cdeltat = 0.1
# number of subdivisions per t=1 unit for bending calculations
deltat = 0.001

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


function fixwidth(xlist, idealwidth)
    initialsep = xlist[end] - xlist[1]
    ratio = idealwidth/initialsep
    newlist = xlist .* ratio
    return newlist
end


#=

Import the initial configuration

=#
# make the b spline basis st the order is 5 and it has C3 continuity
# the second argument determines the range of t
# ALL THE SPLINES HAVE THE SAME BASIS
basis = BSplineBasis(5, 0:5)

# # import the data to make the curve
# initialcurve = CSV.read("filename.csv")

# # make an initial spline to represent the imported curve
# spl = approximate() # okay so I can't get this to work rn so maybe just import the initial control points?
# # alternatively, could use the fitting control points fnc from basicbspline.jl package BUT IT NEEDS A FUNCTION??

# initial control points
# nine for each dimension
xpoints =  [-1.34, -1.13341, -0.512545, -0.125985, -0.690772, -1.24648, -0.846774, -0.223135, 0]
ypoints = [-0.681032, -0.713483, -0.572464, -0.0162971, 0.920295, 1.85729, 2.40882, 2.54411, 2.51039]
zpoints = [0.0, -0.023869, -0.082746, -0.725472, -1.12123, -0.715115, -0.074672, -0.022724,  0.000478]


# scale the x control points so it matches the cell width
# prevents overlapping with the compression fnc
xpoints = fixwidth(xpoints, width/2)


# plot initial curve
# xspline = Spline(basis, xpoints)
# yspline = Spline(basis, ypoints)
# zspline = Spline(basis, zpoints)
# curve = [[xspline(t) yspline(t) zspline(t)] for t in 0:0.01:5]
# plotcurve = plot(Tuple.(curve), xlabel="X", ylabel="Y", zlabel="Z")
# png(plotcurve, "initialcurve.png")




#=

Bending

=#

# function to calculate bending energy
function bending(dcurve, ddcurve, dbasis, ddbasis)
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
function dEbenddPi(lista, tmin, tmax, dimension, knot)
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

        # two in denominator is the 1/2 from B/2 in bending def
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
        f = force(EPS) + 100 * ((dist - (inner_dia + EPS*dia_diff))/(inner_dia + EPS*dia_diff))^2
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
function compression(curve, dcurve, thebasis, dbasis, height, width)
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

    # list = []

    # println("first point in curve")
    # println(curve[1])
    # println("last point in curve")
    # println(curve[end])

    # println("first point in right_10")
    # println(right_10[1])
    # println("last point in right_10")
    # println(right_10[end])

    # find all the contacts
    for i in 1:length(curve)
        normdr1 = norm(dcurve[i])

        # position
        r1x = curve[i][1]
        r1y = curve[i][2]
        r1z = curve[i][3]

        # derivative of the curve along all points
        dr1x = dcurve[i][1]
        dr1y = dcurve[i][2]
        dr1z = dcurve[i][3]

        wtf = []

        # self contacts
        for j in 1:length(curve)
            selfenergy = 0
            t1 = curve[i][4]
            t2 = curve[j][4]
            if abs(t1-t2) > 1 && t1 < t2 #prevents double counting
                normdr2 = norm(dcurve[j])

                # positions
                r2x = curve[j][1]
                r2y = curve[j][2]
                r2z = curve[j][3]

                # derivatives
                dr2x = dcurve[j][1]
                dr2y = dcurve[j][2]
                dr2z = dcurve[j][3]

                # gives the t location for each point on each curve
                t1 = curve[i][4]
                t2 = curve[j][4]

                # distance between the points on each curve
                dist = distance(curve[i], curve[j])

                if dist <= 2.1*ryarn
                    push!(wtf, [dist, t1, t2])
                    contact = floof(dist)

                    # to get potential energy, multiply the energy density by an area and then the coord transform
                    potentialdensity = contact[1]
                    potential = potentialdensity*normdr1*normdr2*cdeltat^2
                    
                    # this forcemag is still a density
                    forcemag = contact[2]

                    # direction of the force
                    rhat = [(r1x-r2x), (r1y-r2y), (r1z-r2z)]
                    rhat = rhat./dist

                    
                    totcompeng += potential
                    selfenergy += potential

                    # for use in computing dE/dP_I
                    # force on left side
                    push!(fordEcomp, [t1,normdr1,dr1x,dr1y,dr1z,t2,normdr2,dr2x,dr2y,dr2z,forcemag,rhat,potentialdensity, 
                    thebasis[i], dbasis[i], thebasis[j], dbasis[j],r1x, r1y, r1z, r2x, r2y, r2z])
                    
                    # DOES INCLUDING THIS FORCE DOUBLECOUNTING??

                    # # force on right side
                    # push!(fordEcomp, [t2,normdr2,dr2x,dr2y,dr2z,t1,normdr1,dr1x,dr1y,dr1z,forcemag,-1*rhat,potentialdensity, 
                    # thebasis[j], dbasis[j], thebasis[i], dbasis[j],r2x, r2y, r2z, r1x, r1y, r1z])
                end
            end
        end

        # contacts with right side
        for j in 1:length(right00)
            if right00[j][4] <= 3.5 #t2<3.5
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
                dist = distance(curve[i], right00[j])

                if dist <= 2.1*ryarn
                    contact = floof(dist)

                    # to get potential energy, multiply the energy density by an area and then the coord transform
                    potentialdensity = contact[1]
                    potential = potentialdensity*normdr1*normdr2*cdeltat^2
                    
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
                    
                    #WOULD THIS CAUSE DOUBLE COUNTING??
                    # # force on right side
                    # push!(fordEcomp, [t2,normdr2,dr2x,dr2y,dr2z,t2,normdr2,dr2x,dr2y,dr2z,forcemag,-1*rhat,potentialdensity, 
                    # thebasis[i], dbasis[i], thebasis[j], dbasis[j],r1x, r1y, r1z, r2x, r2y, r2z])
                end
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
            dist = distance(curve[i], left01[j])

            if dist <= 2.1*ryarn
                contact = floof(dist)

                # to get potential energy, multiply the energy density by an area and then the coord transform
                potentialdensity = contact[1]
                potential = potentialdensity*normdr1*normdr2*cdeltat^2
                
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
                
                # need to fix this here so the thing added to fordEcomp reflects the right thing
                # my current thought is that I need to switch the locations of the contact. so switch t1, 
                # normdr1, dr1x, dr1y, dr1z with there 2 couterparts
                # I already have rhat correctly transformed
                # thebasis and dbasis also need to be switched
                # and then I can do translate on the positions and push them down one cell

                push!(fordEcomp, [t2,normdr2,dr2x,dr2y,dr2z,t1,normdr1,dr1x,dr1y,dr1z,forcemag,minusrhat,potentialdensity, 
                thebasis[j], dbasis[j], thebasis[i], dbasis[i],r2x, r2y-height, r2z, r1x, r1y-height, r1z])
                
                # the 1 down one
                #push!(fordEcomp, [t1,normdr1,dr1x,dr1y,dr1z,t2,normdr2,dr2x,dr2y,dr2z,forcemag,minusrhat,potentialdensity, 
                #thebasis[i], dbasis[i], thebasis[j], dbasis[j],r1x, r1y, r1z, r2x, r2y-2*height, r2z])
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
            dist = distance(curve[i], left02[j])

            if dist <= 2.1*ryarn
                contact = floof(dist)

                # to get potential energy, multiply the energy density by an area and then the coord transform
                potentialdensity = contact[1]
                potential = potentialdensity*normdr1*normdr2*cdeltat^2
                
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

                # need to fix this here so the thing added to fordEcomp reflects the right thing
                # my current thought is that I need to switch the locations of the contact. so switch t1, 
                # normdr1, dr1x, dr1y, dr1z with there 2 couterparts
                # I already have rhat correctly transformed
                # thebasis and dbasis also need to be switched
                # and then I can do translate on the positions and push them down one cell

                push!(fordEcomp, [t2,normdr2,dr2x,dr2y,dr2z,t1,normdr1,dr1x,dr1y,dr1z,forcemag,minusrhat,potentialdensity, 
                thebasis[j], dbasis[j], thebasis[i], dbasis[i],r2x, r2y-2*height, r2z, r1x, r1y-2*height, r1z])

                # the 2 down one
                #push!(fordEcomp, [t1,normdr1,dr1x,dr1y,dr1z,t2,normdr2,dr2x,dr2y,dr2z,forcemag,minusrhat,potentialdensity, 
                #thebasis[i], dbasis[i], thebasis[j], dbasis[j],r1x, r1y, r1z, r2x, r2y-4*height, r2z])
            end
        end

        # contacts with left 1, right side of stitch
        for j in 1:length(right_10)
            if right_10[j][4] >= 1.5 #t2>=1.5
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
                dist = distance(curve[i], right_10[j])

                if dist <= 2.1*ryarn
                    # push!(list,dist)
                    contact = floof(dist)

                    # to get potential energy, multiply the energy density by an area and then the coord transform
                    potentialdensity = contact[1]
                    potential = potentialdensity*normdr1*normdr2*cdeltat^2
                    
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
    end
    # println(length(list))
    # println(minimum(list))
    return totcompeng, fordEcomp
end

# for computing the change in compression energy wrt a control point
# tmin and tmax are the bounds of the parameter affected by the control point
# dimension is x=1,y=2,z=3
# knot indicates which basis function is affected by the control point
function dEcompdPi(list, tmin, tmax, dimension, knot)
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
            dEdP = cdeltat^2*(firstterm + secondterm)

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
function totallength(list, tmax)
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
function dLdPi(list, tmin, tmax, dimension, knot, currentlength)

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
            push!(basislist,[0, autobasis[2], autobasis[3], autobasis[4], autobasis[5], autobasis[6], 0, 0, 0])
        end
    end
    if tmin<=2 && tmax>=3
        for t in 2*units:(3*units-1)
            autobasis = bsplines(basis, deltat*t, Derivative(d))
            push!(basislist,[0, 0, autobasis[3], autobasis[4], autobasis[5], autobasis[6], autobasis[7], 0, 0])
        end
    end
    if tmin<=3 && tmax>=4
        for t in 3*units:(4*units-1)
            autobasis = bsplines(basis, deltat*t, Derivative(d))
            push!(basislist,[0, 0, 0, autobasis[4], autobasis[5], autobasis[6], autobasis[7], autobasis[8], 0])
        end
    end
    if tmin<=4 && tmax>=5
        for t in 4*units:5*units
            autobasis = bsplines(basis, deltat*t, Derivative(d))
            push!(basislist,[0, 0, 0, 0, autobasis[5], autobasis[6], autobasis[7], autobasis[8], autobasis[9]])
        end
    end


    return basislist
end


#=

Stitch Energy

=#

function findbendingenergy(spline1, spline2, spline3)
    # gets values for the curve positions, and 1st 2nd derivatives
    curve = [[spline1(t), spline2(t), spline3(t), t] for t in 0:deltat:5]
    dcurve = [[spline1(t, Derivative(1)), spline2(t, Derivative(1)), spline3(t, Derivative(1))] for t in 0:deltat:5]
    ddcurve = [[spline1(t, Derivative(2)), spline2(t, Derivative(2)), spline3(t, Derivative(2))] for t in 0:deltat:5]

    # gets values of the basis functions and their derivatives at every point
    thebasis = constructbasis5(basis, deltat, 0, 5, 0)
    dbasis = constructbasis5(basis, deltat, 0, 5, 1)
    ddbasis = constructbasis5(basis, deltat, 0, 5, 2)

    #computes bending energy for curve
    bendinge = bending(dcurve, ddcurve, dbasis, ddbasis)
    # println("bending")
    # println(bendinge[1])
    return bendinge
end

function findcontactenergy(spline1, spline2, spline3, height, width)
    # gets values for the curve positions, and 1st 2nd derivatives
    curve = [[spline1(t), spline2(t), spline3(t), t] for t in 0:cdeltat:5]
    dcurve = [[spline1(t, Derivative(1)), spline2(t, Derivative(1)), spline3(t, Derivative(1))] for t in 0:cdeltat:5]
    ddcurve = [[spline1(t, Derivative(2)), spline2(t, Derivative(2)), spline3(t, Derivative(2))] for t in 0:cdeltat:5]

    # gets values of the basis functions and their derivatives at every point
    thebasis = constructbasis5(basis, cdeltat, 0, 5, 0)
    dbasis = constructbasis5(basis, cdeltat, 0, 5, 1)
    ddbasis = constructbasis5(basis, cdeltat, 0, 5, 2)

    # computes compression energy
    contact = compression(curve, dcurve, thebasis, dbasis, height, width)
    # println("contact")
    # println(contact[1])
end

# gives energy of the left half of the stitch
function totalenergy(cpt1, cpt2, cpt3, height, width)
    # contructs the splines
    spline1 = Spline(basis, cpt1)
    spline2 = Spline(basis, cpt2)
    spline3 = Spline(basis, cpt3)

    bendinge = findbendingenergy(spline1, spline2, spline3)
    contact = findcontactenergy(spline1, spline2, spline3, height, width)

    energy = contact[1] + bendinge[1]

    return energy, bendinge[2], contact[2]
end

# gives energy of the left half of the stitch
function totalenergy(cpt1, cpt2, cpt3, height, width, converged)
    # contructs the splines
    spline1 = Spline(basis, cpt1)
    spline2 = Spline(basis, cpt2)
    spline3 = Spline(basis, cpt3)

    bendinge = findbendingenergy(spline1, spline2, spline3)
    contact = findcontactenergy(spline1, spline2, spline3, height, width)

    energy = contact[1] + bendinge[1]

    # prints out the total energy of the entire stitch
    if converged == true
        println("bending")
        println(2*bendinge[1])
        println("contact")
        println(2*contact[1])
        println("total energy")
        println(2*energy)
    end
end

#=

Descent

=#

# HERE IS WHERE THE DESCENT IS DEFINED
# need to do something for the height
function descent(listbend, listcomp, tmin, tmax, dimension, knot, lambda, length)
    benddescent = dEbenddPi(listbend, tmin, tmax, dimension, knot)
    compdescent = dEcompdPi(listcomp, tmin, tmax, dimension, knot)
    lengthstuff = dLdPi(listbend, tmin, tmax, dimension, knot, length)
    lengthascent = lengthstuff[1]
    lengthpenalty = lengthstuff[2]

    # height stuff??
    # need these for all the knots?
    #heightdescent = dEcompdPi(listcomp, deltat, tmin, tmax, 2, knot)


    descent =  compdescent + lambda*lengthascent +  benddescent + penalty*lengthpenalty
    return descent
end


# need to update/fix
function gradientdescent(cpt1, cpt2, cpt3, learn_rate, conv_threshold, max_iter)

    # fix the width of the stitch
    # cpt1 = fixwidth(cpt1, width)

    # gets all the data for the original input curves
    oldenergy = totalenergy(cpt1, cpt2, cpt3, height, width)
    fordEbend = oldenergy[2]
    fordEcomp = oldenergy[3]
    oldenergy = oldenergy[1]

    oglength = totallength(fordEbend,5)

    println("Here is the original energy")
    println(oldenergy)

    # prints the components of the energy
    #printtotalenergy(cpt1, cpt2, cpt3)

    println("here is the og length of curve")
    println(oglength)

    lambda = 0.001

    converged = false
    iterations = 0

    while converged == false

        length = totallength(fordEbend,5)

        # control points for x-direction of curve
        # comment out the descent on the first control pt to fix the width
        cpt101 = cpt1[1] #-learn_rate*descent(fordEbend, fordEcomp, 0, 1, 1, 1, lambda, length)
        cpt102 = cpt1[2] -learn_rate*descent(fordEbend, fordEcomp, 0, 2, 1, 2, lambda, length)
        cpt103 = cpt1[3] -learn_rate*descent(fordEbend, fordEcomp, 0, 3, 1, 3, lambda, length)
        cpt104 = cpt1[4] -learn_rate*descent(fordEbend, fordEcomp, 0, 4, 1, 4, lambda, length)
        cpt105 = cpt1[5] -learn_rate*descent(fordEbend, fordEcomp, 0, 5, 1, 5, lambda, length)
        cpt106 = cpt1[6] -learn_rate*descent(fordEbend, fordEcomp, 1, 5, 1, 6, lambda, length)
        cpt107 = cpt1[7] -learn_rate*descent(fordEbend, fordEcomp, 2, 5, 1, 7, lambda, length)
        cpt108 = cpt1[8] -learn_rate*descent(fordEbend, fordEcomp, 3, 5, 1, 8, lambda, length)
        cpt109 = cpt1[9] #-learn_rate*descent(fordEbend, fordEcomp, 4, 5, 1, 9, lambda, length)
        # last control point of x should be 0 and stay 0

        cpt1 = [cpt101, cpt102, cpt103, cpt104, cpt105, cpt106, cpt107, cpt108, cpt109]

        width = 2*(cpt109 - cpt101)

        # control points for y-direction of curve
        cpt201 = cpt2[1] -learn_rate*descent(fordEbend, fordEcomp, 0, 1, 2, 1, lambda, length)
        cpt202 = cpt2[2] -learn_rate*descent(fordEbend, fordEcomp, 0, 2, 2, 2, lambda, length)
        cpt203 = cpt2[3] -learn_rate*descent(fordEbend, fordEcomp, 0, 3, 2, 3, lambda, length)
        cpt204 = cpt2[4] -learn_rate*descent(fordEbend, fordEcomp, 0, 4, 2, 4, lambda, length)
        cpt205 = cpt2[5] -learn_rate*descent(fordEbend, fordEcomp, 0, 5, 2, 5, lambda, length)
        cpt206 = cpt2[6] -learn_rate*descent(fordEbend, fordEcomp, 1, 5, 2, 6, lambda, length)
        cpt207 = cpt2[7] -learn_rate*descent(fordEbend, fordEcomp, 2, 5, 2, 7, lambda, length)
        cpt208 = cpt2[8] -learn_rate*descent(fordEbend, fordEcomp, 3, 5, 2, 8, lambda, length)
        cpt209 = cpt2[9] -learn_rate*descent(fordEbend, fordEcomp, 4, 5, 2, 9, lambda, length)

        cpt2 = [cpt201, cpt202, cpt203, cpt204, cpt205, cpt206, cpt207, cpt208, cpt209]

        # control points for z-direction of curve
        cpt301 = cpt3[1] #-learn_rate*descent(fordEbend, fordEcomp, deltat, 0, 1, 3, 1, lambda, length)
        cpt302 = cpt3[2] -learn_rate*descent(fordEbend, fordEcomp, 0, 2, 3, 2, lambda, length)
        cpt303 = cpt3[3] -learn_rate*descent(fordEbend, fordEcomp, 0, 3, 3, 3, lambda, length)
        cpt304 = cpt3[4] -learn_rate*descent(fordEbend, fordEcomp, 0, 4, 3, 4, lambda, length)
        cpt305 = cpt3[5] -learn_rate*descent(fordEbend, fordEcomp, 0, 5, 3, 5, lambda, length)
        cpt306 = cpt3[6] -learn_rate*descent(fordEbend, fordEcomp, 1, 5, 3, 6, lambda, length)
        cpt307 = cpt3[7] -learn_rate*descent(fordEbend, fordEcomp, 2, 5, 3, 7, lambda, length)
        cpt308 = cpt3[8] -learn_rate*descent(fordEbend, fordEcomp, 3, 5, 3, 8, lambda, length)
        cpt309 = cpt3[9] -learn_rate*descent(fordEbend, fordEcomp, 4, 5, 3, 9, lambda, length)
        # first control point of z should be 0 and stay 0

        cpt3 = [cpt301, cpt302, cpt303, cpt304, cpt305, cpt306, cpt307, cpt308, cpt309]

        #computes the new energies and the new lists used to determine dE/dP_i
        newenergy = totalenergy(cpt1, cpt2, cpt3, height, width)
        fordEbend = newenergy[2]
        fordEcomp = newenergy[3]
        newenergy = newenergy[1]

        println("here is the updated energy")
        println(newenergy)
        println(iterations)

        #gradient ascent part for the length constraint
        lambda = lambda + (totallength(fordEbend,5) - targetlength)

        if abs(oldenergy - newenergy) <= conv_threshold && iterations > 1000
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
    println(totallength(fordEbend,5))

    println("Here are the new control points")

    println(cpt1)
    println(cpt2)
    println(cpt3)

    println("Here is the energy breakdown:")

    totalenergy(cpt1, cpt2, cpt3, height, width, true)

    return cpt1, cpt2, cpt3
end


#=

Where the code actually runs

=#

grad = gradientdescent(xpoints, ypoints, zpoints, 0.0005, 0.000000005, 100000)


# converged = true

# energy = totalenergy(xpoints, ypoints, zpoints, height, width)
# println(energy[1])

# comptest = compression(curve, dcurve, deltat, thebasis, dbasis, height, width)
# println(comptest[1])
#println(minimum(comptest[2]))

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

curve = [[xspline(t) yspline(t)] for t in 0:0.01:5]
plotcurve = plot(Tuple.(curve), xlabel="X", ylabel="Y")
png(plotcurve, "finalcurvexyplane.png")
fuck = (yspline(5) - yspline(4.99))/(xspline(5)-xspline(4.99))
println(fuck)