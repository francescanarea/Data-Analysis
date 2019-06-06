from sympy import *

print("\n\n\nThis program will create a model for environmental concerns within the animal agriculture industry.")
print("\n\n\nHere are my initial assumptions:")
print("1. Cows that are born to be slaughtered consume water (both directly and indirect, such as through eating wheat)")
print("2. With not enouhg water in the agriculture system, the number of cows enslaved for purposes of human consumption will decrease (these cows will either die due to lack of nutrition based on low amounts of food supply or water to keep them alive and, also, longerterm slaughter-houses will not force as much procreation based on their lack of resources so less enslaved cows will be born).")
print("3.With not enough cows enslaved and to be slaughtered, then the water supply within the agriculture system will increase.\n\n\nSummary:\n1.Cows consume water.\n2.Less water = less cows.\n3.Less cows = more water")

print("\n\n\nWe will let Xk = [Ck, Wk], where\nCk = num of cows in millions to be slaughtered at year t=k")
print("Wk = num olympic swimming pools of water in agricultural system")


print("\n\n\nHere's some initial data.")

print("\nCow data:")
#cows in millions
cows_2014 = 300.07
cows_2004 = 279.91

print("Cows slaughtered in 2004: " + str(cows_2004))
print("Cows slaughtered in 2014: " + str(cows_2014))

change_cows = cows_2014 - cows_2004
increase_cows = (change_cows/cows_2004)/10

print("2004-2014 increase of cows slaughtered: " + str(increase_cows))

print("\nWater data:")
#water measured as the number of olympic swimming pools
gal_per_pool = 660253.09 #660,253.09 galloons/olympic swimming pool

kg_per_mil_cows_2004 = 202.6*1000000
kg_per_mil_cows_2014 = 209.6*1000000

liters_per_kg = 15415

liters_per_gal = 3.78541  #x liters / 15415 = num gal

print("Each kg of cow meat takes " + str(liters_per_kg) + " liters of water to produce.")
print("In 2004, each cow produced 202.6 kg of meat, so 1 million cows produced " + str(kg_per_mil_cows_2004) + " kg of meat.")
print("In 2014, each cow produced 209.6 kg o fmeat, so 1 million cows produced " + str(kg_per_mil_cows_2014) + " kg of meat.") 

#2004-- 1. find num liters, 2. convert to swimming pools
total_kg_2004 = cows_2004*kg_per_mil_cows_2004

total_liters_2004 = total_kg_2004 * liters_per_kg

total_gal_2004 = total_liters_2004 / liters_per_gal

pools_2004 = total_gal_2004 / gal_per_pool

#2014
total_kg_2014 = cows_2014*kg_per_mil_cows_2014

total_liters_2014 = total_kg_2014 * liters_per_kg

total_gal_2014 = total_liters_2014 / liters_per_gal

pools_2014 = total_gal_2014 / gal_per_pool

change_pools = pools_2014 - pools_2004
increase_pools = (change_pools / pools_2004)/10

print("\nIn 2004, the totals for cow meat production were: " + str(int(total_kg_2004)) + " kg of meat.")
print("Thus, the total liters of water used for cow meat production in 2004 was: " + str(int(total_liters_2004)) + ", which translates to " + str(pools_2004) + " olympic swimming pools filled with water.")

print("\nIn 2014, the totals for cow meat production were: " + str(total_kg_2014))
print("So, in 2014, the total liters of water used for cow meat production were: " + str(total_liters_2014) + ", which translates to " + str(pools_2014) + " olympic swimming pools.")

print("So the increase of water consumption from 2004-2014 was: " +str(increase_pools) + " swimming pools.")

print("\n\n\nFrom here, we can begin creating our equations: \n\n\n")

print("Ck+1= " + str(increase_cows) + "Ck + " + str(increase_pools) + "Wk")

print("Wk+1 = " + str(increase_cows*-1) + "Ck + " + str(1+increase_pools) + "Wk")

A = Matrix([[increase_cows, increase_pools], [increase_cows*-1, 1+increase_pools]])

print("Our matrix to model these equations is: " + str(A))

print("\n\n\nThe first step is to compute the eigenvalues of this matrix.  We begin by solving A-(lambda)I.")
x = symbols('x')
I = eye(2)
lambda_i = x*I
print("\nSolved (lambda)*I: " + str(lambda_i))

P = A - lambda_i
print("\nSolved A-(lambda)*I:" + str(P))

eq1 = expand((P.row(0).col(0))*(P.row(1).col(1)))
#print(eq1)

poly = Poly(eq1[0], x).all_coeffs()
#print("\nPolynomial: " + str(poly))

#where poly[0] = coeff of x**2
a = poly[0]
b = poly[1]
c = poly[2]

print("\n\nSolve det(A-(lambda)I): " + str(poly))
print("\nNext, we need to use the quadratic equation to find the roots of this polynomial.")
print("\nMy a= " + str(a) + ", b = " + str(b) + " and c = " + str(c))
lambda_1 = (-1*b +sqrt(b**2 - 4*a*c))/(2*a)
lambda_2 = -1*b - sqrt(b**2 - 4*a*c)/(2*a)

print("Using quadratic equation to find eigenvalues...\n\n\n")
print("Success!\n\n\n")
print("\nLambda 1 = " + str(lambda_1) + "\nLambda 2 = " + str(lambda_2))

#EIGENVALUES
a1 = (A - lambda_1*I) 
a2 = (A - lambda_2*I)

print("\n\n\nFinding eigenvectors...")
print("\n\nSolving homogeneous equation for: \n" + "(A-" + str(lambda_1) + "*I)x = 0\n(A-" + str(lambda_2) + "*I)x = 0") 

# array of tuples, where eigenvectors[0] = (eigenvalue, num, sympy matrix with eigenvector)
eigenvectors = A.eigenvects()
ev1 = eigenvectors[0][2] #retrieving the first tuple and 3rd value of the tuple which holds the sympy eigenvector matrix
ev2 = eigenvectors[1][2] 

ev1 = ev1[0].tolist() #converting to a regular array so can reference each value of vector easily
ev2 = ev2[0].tolist()

print("\n\n\nEigenvector #1: " + str(ev1) + "\nEigenvector #2: " + str(ev2))

print("\n\n\nFinding dominant eigenvector...")

#finding dominant eigenvector by comparing values of each eigenvector
dominant_eig = ev1 if ((ev1[0] >= ev2[0]) and (ev1[1] >= ev2[1])) else ev2
print("Dominant eigenvector: " + str(dominant_eig) + "\n\n")

f= open("eigenvectors.txt","w+")

f.write(str(ev1[0][0]))
f.write("\n"+str(ev1[1][0]))
f.write("\n"+str(ev2[0][0]))
f.write("\n"+str(ev2[1][0]))

f.close()

