import numpy as np


def solve(f,a: float, b: float, h: float,k: float,conditions_lim: list):
    if a==0 or b==0 or h==0 or k==0:
        return None

    if not conditions_lim:
        return None

    n,m=1+a/h,1+b/k
    dim=(n-2)*(m-2)

    if not (n.is_integer() or n.is_integer()):
        return None

    A=_get_matrix(f,n,m,h,k)
    b=_get_b(f,a,b,dim,conditions_lim)
    
    return None


def _get_matrix(n,m,h,k):
    taille=(n-2)*(m-2)
    A=np.zeros((taille,taille))
    line=0
    for i in range(2,m):
        for j in range(2,n):
            A[line,(i-2)*(n-2)+j-2]=-2*(1/h**2+1/k**2) #coeff de u_ij
            if i>2:
                A[line,(i-3)*(n-2)+j-2]=1/h**2 # coeff de u_{i-1}_j
            if i<m-1:
                A[line,(i-1)*(n-2)+j-2]=1/h**2      #coeff  u_{i+1}_j
            if j>2:
                A[line, (i - 2) * (n - 2) + j - 3] =1/k**2 # coeff u_i_{j-1}
            if j<n-1:
                A[line, (i - 2) * (n - 2) + j - 1] =1/k**2 # coeff u_i_{j+1}

            line+=1

    return A


def _get_b(f,n,m,h,k,conditions_lim):
    dim=(n-2)*(m-2)

    def _get_lim_elm(i,j):
        """
        Prends les coordonees d'un point du bord de la maille et retourne sa valeur (conditions aux limites)
        :param i: ligne i de la maille
        :param j: ligne j de la maille
        :return:
        """

        if i==1:
            return conditions_lim[j-1]
        if i==m:
            return conditions_lim[m+2*n-j-2]
        if j==1:
            return conditions_lim[2*n+2*m-i-3]
        if j==n:
            return conditions_lim[n+i-2]

    bo=np.zeros(dim)
    line=0
    for i in range(2,m):
        for j in range(2,n):
            bo[line] = f(h*(i-1), k*(j-1))
            line+=1
    i_2,i_m_1=2,m-1
    for j in range(2,n):
        # line=(i-2)*(n-2)+j-2 # on est sur la line-ieme equation / ligne dans A avec i,j l'element du maillage i,j>0
        bo[(i_2 - 2) * (n-2) + j - 2] -=_get_lim_elm(1,j)/h**2
        bo[(i_m_1 - 2) * (n-2) + j - 2] -= _get_lim_elm(m,j)/h**2

    j_1,j_n_1=2,n-1
    for i in range(2,m):
        bo[(i - 2) * (n-2) + j_1 - 2] -=_get_lim_elm(i,1)/k**2
        bo[(i - 2) * (n-2) + j_n_1 - 2] -=_get_lim_elm(i,n)/k**2

    return bo




#print(_get_matrix(5,5,1,1))
print(_get_b(lambda x,y:x**2,5,5,1,1,[1]*16))