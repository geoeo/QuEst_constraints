from sympy import *
from sympy.combinatorics import Permutation


x,y,z,w = symbols('x y z w', real=True)

flx1, fly1 = symbols('flx1 fly1', real=True)
flx2, fly2 = symbols('flx2 fly2', real=True)
flx3, fly3 = symbols('flx3 fly3', real=True)
flx4, fly4 = symbols('flx4 fly4', real=True)
flx5, fly5 = symbols('flx5 fly5', real=True)

frx1, fry1 = symbols('frx1 fry1', real=True)
frx2, fry2 = symbols('frx2 fry2', real=True)
frx3, fry3 = symbols('frx3 fry3', real=True)
frx4, fry4 = symbols('frx4 fry4', real=True)
frx5, fry5 = symbols('frx5 fry5', real=True)

tx,ty,tz = symbols('tx ty tz', real=True)
u1,u2,u3,u4,u5 = symbols('u1 u2 u3 u4 u5', real=True)
v1,v2,v3,v4,v5 = symbols('v1 v2 v3 v4 v5', real=True)

feat_left_1 = Matrix([flx1,fly1,-1])
feat_left_2 = Matrix([flx2,fly2,-1])
feat_left_3 = Matrix([flx3,fly3,-1])
feat_left_4 = Matrix([flx4,fly4,-1])
feat_left_5 = Matrix([flx5,fly5,-1])

feat_right_1 = Matrix([frx1,fry1,-1])
feat_right_2 = Matrix([frx2,fry2,-1])
feat_right_3 = Matrix([frx3,fry3,-1])
feat_right_4 = Matrix([frx4,fry4,-1])
feat_right_5 = Matrix([frx5,fry5,-1])

R = Matrix([[Pow(w,2)+Pow(x,2)-Pow(y,2)-Pow(z,2),2*(x*y-w*z),2*(x*z+w*y)],
            [2*(x*y+w*z),Pow(w,2)-Pow(x,2)+Pow(y,2)-Pow(z,2),2*(y*z-w*x)],
            [2*(x*z-w*y),2*(y*z+w*x),Pow(w,2)-Pow(x,2)-Pow(y,2)+Pow(z,2)]])

T = Matrix([[tx],[ty],[tz]])

e1_left = u1*R*feat_left_1+T
e1_right = v1*feat_right_1

e2_left = u2*R*feat_left_2+T
e2_right = v2*feat_right_2

e3_left = u3*R*feat_left_3+T
e3_right = v3*feat_right_3

e4_left = u4*R*feat_left_4+T
e4_right = v4*feat_right_4

e5_left = u5*R*feat_left_5+T
e5_right = v5*feat_right_5

M = Matrix.zeros(6,6)
M[0:3,0] = R*feat_left_1
M[0:3,1] = feat_right_1
M[0:3,2] = R*feat_left_2
M[0:3,3] = feat_right_2

M[3:6,0] = R*feat_left_1
M[3:6,1] = feat_right_1
M[3:6,4] = R*feat_left_3
M[3:6,5] = feat_right_3

print("Calc Determinant...")
#det_m = det(M)

(L,U,perm) = M.LUdecomposition()
P = eye(M.rows).permuteFwd(perm)
P_inv = P.inv()

det_P_inv = det(P_inv)
det_L = det(L)
det_U = det(U)

det_M = det_P_inv*det_L*det_U

print(det_M)

