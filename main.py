import math
import numpy as np
import cmath
from sklearn.decomposition import PCA, FastICA
from sklearn.preprocessing import StandardScaler

def new_bar(matriz):
    #cria barras. [tipo de barra(0), capacitância(1), tensão(2), angulo(3), potência ativa(4), potência reativa(5), existe carga?(6),retância positiva(7), tensão base(8), strafo(9)]
    #tensão em p.u., ângulo em radianos, potência ativa e reativa em p.u.

    bar = [0,0,0,0,0,0,0,0,0,0]

    bartype = int(input("Insira o tipo da barra: "))
    carga = int(input("Existe carga? (1 = sim, 0 = não): "))
    gerador = float(input("Insira o valor de reatância positiva do gerador (p.u.): "))
    c = float(input("Insira a potência reativa da capacitância da barra (p.u.): "))
    vb = float(input("Insira o valor da tensão base da barra (kV): "))
    t = int(input("Existe trafo? (1 = sim, 0 = não): "))
    if(t == 1):
        st = float(input("Insira o valor Strafo (p.u.): "))
    else:
        st = 0

    if bartype == 1:

        v = float(input("Insira a tensão (p.u.): "))
        dg = float(input("Insira o angulo (rad): "))
        bar = [bartype,c,v,dg,0,0,carga, gerador, vb, st]

    else:
        if bartype == 2:
            v = float(input("Insira a tensão (p.u.): "))
            p = float(input("Insira a potência ativa gerada (pu): "))

            if(carga == 0):
                bar = [bartype,c,v,0,p,0,carga, gerador, vb, st]
            else:
                bar = [bartype,c,v,0,-p,0,carga, gerador, vb, st]
        else:
            if bartype == 3:
                p = float(input("Insira a potência ativa liberada (pu): "))
                q = float(input("Insira a potência reativa liberada (pu): "))

                if(carga == 0):
                    bar = [bartype,c,1,0,p,q,carga, gerador, vb, st]
                else:
                    bar = [bartype,c,1,0,-p,-q,carga, gerador, vb, st]

    matriz.append(bar)

def new_connect(y, frombar, tobar,bars, connections):
    #cria linhas
    r = float(input("Insira a resistência (p.u.): "))
    xl = float(input("Insira a reatância indutiva (p.u.): "))
    b = float(input("Insira a susceptância (p.u.): "))

    addconnect = [frombar, tobar, r, xl, b]
    connections.append(addconnect)

    z = r + (xl*1j) 

    y[frombar-1][tobar-1] += -1/z
    y[tobar-1][frombar-1] += -1/z

    y[frombar-1][frombar-1] += 1/z + (b*1j)/2
    y[tobar-1][tobar-1] += 1/z + (b*1j)/2   

def soma_shunt(y, bars):
    #soma elemento shunt das barras
    for i in range(len(bars)):
        if(bars[i][1] != 0):
            y[i][i] += (bars[i][1]*1j)


def new_matriz(tam1, tam2):
    #criação genérica de matriz tam1xtam2
    matriz = []
    for l in range(tam1):
        linha = []
        for c in range(tam2):
            elemento = 0
            linha.append(elemento)
        matriz.append(linha)
    return matriz

def convergence(x, e):
    #teste de convergência
    for i in range(len(x)):
        if(abs(x[i]) > e):
            return False
    
    return True

def bars_qnt(bars):
    #quantidade de PQ e PV (PQ é qnt[0], PV é qnt[1])
    qnt = [0, 0]
    for i in range(len(bars)):
            if bars[i][0] == 2:
                qnt[0] = qnt[0] + 1
            
            if bars[i][0] == 3:
                qnt[1] = qnt[1] + 1      

    return qnt

def create_x(bars, Y, qnt):
    #matriz de delta P e delta Q
    x = new_matriz(qnt[0] + (2 * qnt[1]), 1)

    i = 0

    for i1 in range(len(bars)):
        if bars[i1][0] != 1:
            dP = 0
            for i2 in range(len(Y[0])):
                if Y[i1][i2] != 0:
                    dP += bars[i2][2] * (Y[i1][i2].real * math.cos(bars[i1][3] - bars[i2][3]) + Y[i1][i2].imag * math.sin(bars[i1][3] - bars[i2][3]))
            

            x[i] = bars[i1][4] - (bars[i1][2] * dP)
            i += 1

    for i1 in range(len(bars)):
        if bars[i1][0] == 3:
            dQ = 0
            for i2 in range(len(Y[0])):
                if Y[i1][i2] != 0:
                    dQ += bars[i2][2] * (Y[i1][i2].real * math.sin(bars[i1][3] - bars[i2][3]) - Y[i1][i2].imag * math.cos(bars[i1][3] - bars[i2][3]))

            x[i] = bars[i1][5] - (bars[i1][2] * dQ)
            i +=  1

    return x

def create_h(bars, Y):
    #Cria submatriz h
    H = new_matriz(len(bars) - 1, len(bars) - 1)

    b = 0
    b2 = 0
    

    for i in range(len(H)):
        select = 0
        j = 0
        while(j<len(H[0])):
            valor = 0
            while(select == 0):
                if bars[b][0] != 1:
                    select = 1
                else:
                    b+=1

            for b2 in range(len(Y)):
                valor = 0
                if bars[b2][0] != 1:
                    if b == b2:
                        for n in range(len(Y)):
                            if Y[b][n] != 0 and b != n:
                                valor += bars[n][2] * (-Y[b][n].real * math.sin(bars[b][3] - bars[n][3]) + Y[b][n].imag * math.cos(bars[b][3] - bars[n][3]))
                        H[i][j] = -(valor * bars[b][2])

                    else:
                            valor = bars[b][2]*bars[b2][2] * (Y[b][b2].real * math.sin(bars[b][3] - bars[b2][3]) - Y[b][b2].imag * math.cos(bars[b][3] - bars[b2][3]))
                            H[i][j] = -(valor)
                    j+=1      

            
        b+=1
                
    return H

def create_n(bars, Y):
    #Cria submatriz n
    qnt = bars_qnt(bars)
    N = new_matriz(len(bars) - 1, qnt[1])

    b = 0
    b2 = 0

    for i in range(len(N)):
        select = 0
        j = 0
        while j < len(N[0]):
            valor = 0
            while(select == 0):
                if bars[b][0] != 1:
                    select = 1
                else:
                    b+=1
            
            for b2 in range(len(Y)):
                valor = 0
                if bars[b2][0] == 3:
                    if b == b2:
                        for n in range(len(Y)):
                            if Y[b][n] != 0 and b != n:
                                valor += bars[n][2] * (Y[b][n].real * math.cos(bars[b][3] - bars[n][3]) + Y[b][n].imag * math.sin(bars[b][3] - bars[n][3]))

                        if valor != 0:
                            N[i][j] = -(valor + 2*bars[b][2]*Y[b][b].real)

                    else:
                            valor = bars[b][2] * (Y[b][b2].real * math.cos(bars[b][3] - bars[b2][3]) + Y[b][b2].imag * math.sin(bars[b][3] - bars[b2][3]))
                            N[i][j] = -(valor) 
                    j+=1 
        b+=1
            
    return N

def create_m(bars, Y):
    #Cria submatriz m
    qnt = bars_qnt(bars)
    M = new_matriz(qnt[1], len(bars) - 1)

    b = 0
    b2 = 0

    for i in range(len(M)):
        select = 0
        j = 0
        while j <len(M[0]):
            valor = 0
            while(select == 0):
                if bars[b][0] == 3:
                    select = 1
                else:
                    b+=1

            for b2 in range(len(Y)):
                valor = 0
                if bars[b2][0] != 1:
                    if b == b2:
                        for n in range(len(Y)):
                            if Y[b][n] != 0 and b!=n:
                                        valor += bars[n][2] * (Y[b][n].real * math.cos(bars[b][3] - bars[n][3]) + Y[b][n].imag * math.sin(bars[b][3] - bars[n][3]))
                        M[i][j] = -(valor * bars[b][2])
                    else:

                        valor = -bars[b][2] * bars[b2][2] * (Y[b][b2].real * math.cos(bars[b][3] - bars[b2][3]) + Y[b][b2].imag * math.sin(bars[b][3] - bars[b2][3]))
                        M[i][j] = -(valor)
                    j+=1
        b+=1
    return M

def create_l(bars, Y):
    #Cria submatriz l
    qnt = bars_qnt(bars)
    L = new_matriz(qnt[1], qnt[1])

    b = 0
    b2 = 0

    for i in range(len(L)):
        select = 0
        j = 0
        while j < len(L[0]):
            valor = 0
            while(select == 0):
                if bars[b][0] == 3:
                    select = 1
                else:
                    b+=1  
            
            for b2 in range(len(Y)):
                valor = 0
                if bars[b2][0] == 3:
                    if b == b2:
                        for n in range(len(Y)):
                            if Y[b][n] != 0 and b!=n:
                                        valor += bars[n][2] * (Y[b][n].real * math.sin(bars[b][3] - bars[n][3]) - Y[b][n].imag * math.cos(bars[b][3] - bars[n][3]))
                        
                        if valor!=0:
                            L[i][j] = -(valor + -2*bars[b][2]*Y[b][b].imag)

                    else:
                        valor = bars[b][2] * (Y[b][b2].real * math.sin(bars[b][3] - bars[b2][3]) - Y[b][b2].imag * math.cos(bars[b][3] - bars[b2][3]))
                        L[i][j] = -(valor)
                    j+=1   
        b+=1  

            
                
    return L

def create_jacob(bars, Y, qnt):
    #Une as submatrizes
    H = create_h(bars,Y)
    N = create_n(bars,Y)
    M = create_m(bars,Y)
    L = create_l(bars,Y)

    jacob1 = new_matriz(len(H), len(H[0]) + len(N[0]))

    i1 = 0
    j1 = 0

    for i in range(len(jacob1)):
        j1=0
        for j in range(len(jacob1[0])):
            if(j < len(H[0])):
                jacob1[i][j] = H[i][j]
            else:
                jacob1[i][j] = N[i1][j1]
                j1+=1
        i1+=1
    jacob2 = new_matriz(len(M), len(M[0]) + len(L[0]))

    i1 = 0
    j1 = 0

    for i in range(len(jacob2)):
        j1=0
        for j in range(len(jacob2[0])):
            if(j < len(M[0])):
                jacob2[i][j] = M[i][j]
            else:
                jacob2[i][j] = L[i1][j1]
                j1+=1
        i1+=1

    result = np.vstack((jacob1, jacob2))

    return result
            
def calc_pow(bars,y, qnt):
    #calculo da potência final
    v1 = v2 = 0
    for i in range(len(bars)):
        if bars[i][0] == 1:
            for b in range(len(y)):
                if y[i][b] != 0:
                    v1 += bars[b][2] * (y[i][b].real * math.cos(bars[i][3] - bars[b][3]) + y[i][b].imag * math.sin(bars[i][3] - bars[b][3]))
                    v2 += bars[b][2] * (y[i][b].real * math.sin(bars[i][3] - bars[b][3]) - y[i][b].imag * math.cos(bars[i][3] - bars[b][3]))
            bars[i][4] = v1*bars[i][2]
            bars[i][5] = v2*bars[i][2]
    
        v1 = v2 = 0

        if bars[i][0] == 2:
            for b in range(len(y)):
                if y[i][b] != 0:
                    v2 += bars[b][2] * (y[i][b].real * math.sin(bars[i][3] - bars[b][3]) - y[i][b].imag * math.cos(bars[i][3] - bars[b][3]))
            bars[i][5] = v2*bars[i][2]
                
                
def bars_print(bars):
    #printa as barras em sequência
    for i in range(len(bars)):
        print(bars[i])


def NewtonRhapson(y,bars,qnt,e):
    #Processo do Newton Rhapson de fato
    test = 0
    while True:
        x = create_x(bars, y, barsqnt)
        jacob = np.linalg.inv(np.negative(create_jacob(bars, y,qnt)))

        #print(x)

        result = np.matmul(jacob,x)

        i = 0
        i1 = 0
        #print(create_jacob(bars, y),"\n")
        #bars_print(bars)
        #print(result,"\n")

        while i < qnt[0] + qnt[1]:
            if bars[i1][0] != 1:
                bars[i1][3] = bars[i1][3] + result[i]
                i1+=1
                i+=1
            else:
                i1+=1
        
        i1 = 0

        while  i < (qnt[0] + qnt[1]*2):
            if bars[i1][0] == 3:
                bars[i1][2] = bars[i1][2] + result[i]
                i1+=1
                i+=1
            else:
                i1+=1

        #print(result)
        x = create_x(bars, y, barsqnt)
        #print("test")

        if convergence(x,e) == True:
            calc_pow(bars,y, barsqnt)
            bars_print(bars)
            break

def exist_in(x,vet):
    #verifica se o elemento x existe em um vetor
    flag = False
    count = int
    count = 0
    while count < len(vet):
        if vet[count][0] == x:
            flag = True
        count += 1

    return flag

def get_corrente(barqnts,Sb):
    #recebe as correntes harmôncias
    vet = []
    bar = 1
    while(bar != -1):

        
        bar = int(input("Insira a barra da inserção (-1 para parar): "))
        if(bar != -1):

            h = int(input("Insira a ordem harmônica: "))
            i = float((input("Insira o valor da corrente: ")))/( (Sb*1000)/(bars[bar-1][8]*np.sqrt(3)) )
            f = np.radians((float(input("Insira o fasor da corrente: "))))

            a = ( i * np.cos(f))
            b = ( i * np.sin(f))
            if(exist_in(h,vet)):
                for count in range(len(vet)):
                    if vet[count][0] == h:
                        for count1 in range(len(vet[count][1])):
                            if count1 == bar-1:
                                vet[count][1][count1] = [complex(a,b)]
            else:
                veti = new_matriz(barqnts, 1)
                addvet = [0,0]
                veti[bar-1] = [complex(a,b)]
                addvet[0] = h
                addvet[1] = veti
                #print(addvet,"\n")
                vet.append(addvet)

        #print(vet,"\n")
        
    return vet

def correct_y(connections, barsqnt, h,bars):
    y = new_matriz(barsqnt, barsqnt)
    for count in range(len(connections)):

        frombar = connections[count][0]
        tobar = connections[count][1]
        r = connections[count][2]

        xl = connections[count][3]*h
        b = connections[count][4]*h

        z = r + (xl*1j) 


        y[frombar-1][tobar-1] = -1/z
        y[tobar-1][frombar-1] = -1/z

        y[frombar-1][frombar-1] += 1/z + (b*1j)/2
        y[tobar-1][tobar-1] += 1/z + (b*1j)/2


    for i in range(len(bars)):
        if(bars[i][7] > 0):
            y[i][i] += 1/(bars[i][7]*h*1j)
        
        if(bars[i][6] == 1):
            
            r = (bars[i][2]*bars[i][2])/(-bars[i][4])
            xl = (bars[i][2]*bars[i][2])/(bars[i][5]*1j)
            y[i][i] += 1/r + 1/(xl*h)

        if(bars[i][1] > 0):
            y[i][i] += (bars[i][1]*1j*h)
            #y[i][i] += -(bars[i][2]*bars[i][2])/(bars[i][1]*1j*h)
            #y[i][i] += -((1j*bars[i][1]*h)/(bars[i][2]*bars[i][2]))
    
    #if(h == 3):
        #print(y)
    
    return y


def harmonic_calc(currents, connections, barsqnt,bars,vh):
    DTT = new_matriz(barsqnt, 1)
    DIT = new_matriz(barsqnt, len(currents))
    for count in range(len(currents)):
        h = currents[count][0]

        Yh = np.linalg.inv(correct_y(connections,barsqnt,h,bars))
        #print(currents[count][1])
        #print("\n")
        #print(Yh)

        v = np.matmul(Yh,currents[count][1])

        for l in range(len(DIT)):
            DIT[l][count] = v[l][0]

        for i in range(len(DTT)):
            DTT[i][0] += abs(v[i][0])*abs(v[i][0])

    for i in range(len(DTT)):
        DTT[i][0] = (np.sqrt(DTT[i][0])*100)
    vh.append(v)

    print("MEU VPAC: ",v[1])
    print("DIT: \n",DIT,"\n")
    print("DTT: \n",DTT,"\n")
    return DTT

def calc_impedancias(PAC,bars,h,connections):
    impedancias = new_matriz(1,len(bars)-1)[0]
    fill = 0
    for i in range(len(bars)):
        if(bars[i][0] == 1):
            for c in range(len(connections)):
                if(connections[c][0] == (i+1) or connections[c][1] == (i+1)):
                    #print(c)
                    Rcon = connections[c][2]
                    Xcon = connections[c][3]*h
                    print("Rcon, Xcon:",Rcon,Xcon)
            #Lcon = (Xcon)/(2*math.pi*60)
            Zcon = complex(Rcon,Xcon)*((13800*13800)/10000000)
            impedancias[fill] = Zcon
            #print("Zcon = ",Zcon)
            fill += 1
        else:
            if(i != PAC):
                #Zbase = (v*v)/(bars[i][9]*1000000)
                #print("Zb[",i,"]",Zbase)
                #Rind = ((bars[i][2]*bars[i][8]*1000)*(bars[i][2]*bars[i][8]*1000))/(bars[i][4]*1000000)
                Rind = ((bars[i][2]*13800)*(bars[i][2]*13800))/(bars[i][4]*10000000)
                Rind = Rind*-1
                #print("R[",i,"]",Rind)
                Lind = ((bars[i][2]*13800)*(bars[i][2]*13800)*h)/(bars[i][5]*10000000)
                Lind = Lind*-1
                #print("L[",i,"]",Lind)
                Cind = ((bars[i][2]*13800)*(bars[i][2]*13800))/(bars[i][1]*10000000*h)
                #print("C[",i,"]",Cind)
                #Rind_ = Rind*((bars[PAC][8]/bars[i][8])*(bars[PAC][8]/bars[i][8]))
                #Lind_ = Lind*((bars[PAC][8]/bars[i][8])*(bars[PAC][8]/bars[i][8]))
                #Cind_ = Cind*((bars[i][8]/bars[PAC][8])*(bars[i][8]/bars[PAC][8]))
                for c in range(len(connections)):
                    if(connections[c][0] == (i+1) or connections[c][1] == (i+1)):
                        Rtrafo = connections[c][2]
                        Xtrafo = connections[c][3]*h

                Ztrafo = complex(Rtrafo,Xtrafo)*((13800*13800)/10000000)
                #print("RL_par[",i,"]",RL_par)
                Zind = Ztrafo + 1/(1/Rind + 1/(1j*Lind) + 1/(-1j*Cind))
                #print("Zc[",i,"]",Zc)

                #print("RL_C_par[",i,"]",RL_C_par)
                #print("Z_ind[",i,"]",Z_ind)
                impedancias[fill] = Zind
                if(i == 2):
                    print("Rind: ",Rind)
                    print("Lind: ",Lind)
                    print("Cind: ",Cind)
                    print("Zind[",i,"]",Zind)
                    print("V = ",(bars[i][2]*13800))
                fill += 1

    return impedancias

def calc_impedancias_trans(PAC,c,connections,h):
    impedancias = [0,0]
    for count in range (len(connections)):
        if ((connections[count][0] == c+1) or (connections[count][1] == c+1)):
            impedancias[0] += complex(connections[count][2], connections[count][3]*h)
        else:
            if ((connections[count][0] != c+1) and (connections[count][1] != c+1)):
                impedancias[1] += complex(connections[count][2],connections[count][3]*h)
    return impedancias


def projection_complex(z, w):
    # Cálculo do produto escalar z * w
    dot_product = z.real * w.real + z.imag * w.imag
    
    # Cálculo do módulo de w
    modulus_w = abs(w)
    
    # Cálculo da projeção vetorial
    projection = (dot_product / modulus_w**2) * w
    
    return projection

def defasar_30(z):
    return cmath.rect(abs(z),cmath.phase(z) + (30*(cmath.pi/180)))

def compartilha(y,bars, connections, currents, vh):
    c = int(input("Insira a barra consumidora: "))-1
    for count in range(len(y[c])):
        if y[c][count] != 0 and count != c:
            PAC = count
    zsh = 0
    zch = 0
    impedancias = calc_impedancias(PAC,bars,5,connections)

    #improviso
    impedancias = (cmath.rect(5.9478,(89.5998*(math.pi/180))),cmath.rect(62.5424,(-4.3393*(math.pi/180))),cmath.rect(51.2691,(-4.9947*(math.pi/180))),cmath.rect(83.2596,(-21.6006*(math.pi/180))))
    print(impedancias)
    flag = 0
    # i = 0 flag = 0: zsh += 1/impedancias[0] flag = 1
    # i = 1 flag = 1: i = 2
    # i = 2 flag = 1: zch = impedancias[1] flag = 2
    # i = 3 flag = 2: zsh += 1/impedancias[2] flag = 3
    # i = 4 flag = 3: zsh += 1/impedancias[3] flag = 4
    for i in range(len(bars)):
        #print("i, flag: ",i,flag)
        if(i != PAC and i != c):
            print("Flag, impedância",flag,impedancias[flag])
            zsh += 1/impedancias[flag]
            flag += 1
        else:
            if (i == c):
                zch = impedancias[flag]
                flag += 1
    
    print("Zch: ",abs(zch))
    zsh = 1/zsh
    #zsh = 1.155864 + 5.945482j
    #zsh = zsh/19.044
    #zch = zch/19.044
    print("Zsh: ",abs(zsh))
    print("Zch 19.044: ",abs(zch)*19.044)
    #zch = complex(62.3631,-4.7322)
    #zsh = complex(1.1657,5.9440)
    #trocar 0 para h no final
    print("vh[0][PAC]:",vh[0][PAC])
    print("vh[0][c]:",vh[0][c])

    vc = complex(vh[0][c].real*-1,vh[0][c].imag)
    vpac =  0.00501242+0.032241j
    #ipac = ((vpac - ((vh[0][c]*13800))/math.sqrt(3))/zsh)/math.sqrt(3)
    ipac = ((vpac - ((vc)))/calc_impedancias_trans(PAC, c ,connections, 5)[0])
    ipac = ipac*418.3698
    #ipac = 0.0564415859975+0.0350785718963j
    vpac = (vpac * 13800)

    #print("Defasar 30: ", abs(defasar_30(vpac)), cmath.phase(defasar_30(vpac)))
    ipac = defasar_30(ipac)
    vpac = defasar_30(vpac)
    vpac = vpac/math.sqrt(3)
    vpac = (cmath.rect(264.2867,(128.3914*(math.pi/180))))
    #ipac = -17.52 - 0.434j
    #vpac = -93.853 + 242.431j
    #ipac = ipac/math.sqrt(3)
    #ipac = (ipac * ((10000000)/(13800*np.sqrt(3))))/math.sqrt(3)
    ish = (((vpac)/zsh) + ipac)
    ich = (((vpac)/zch) - ipac)
    #ish = defasar_30(ish)
    #ich = defasar_30(ich)
    #ish = 13.482 + 23.453j
    #ich = 11.995 + 12.823j
    #zsh = (zsh*19.044)
    #zch = (zch*19.044)
    #ish = (ish*418.3698)
    #ich = (ich*418.3698)

    vspach = ((zsh*zch)/(zsh+zch))*ish

    vcpach = ((zsh*zch)/(zsh+zch))*ich
    print("ish: ",abs(ish),cmath.phase(ish)*(180/math.pi) ," -> ",ish)
    print("ich: ",abs(ich),cmath.phase(ich)*(180/math.pi)," -> ",ich)
    print("ipac: ",abs(ipac),cmath.phase(ipac)*(180/math.pi)," -> ",ipac)
    #print("zch: ",abs(zch)*((380*380)/10000000)," -> ",zch)
    print("zsh: ",abs(zsh),cmath.phase(zsh)*(180/math.pi)," -> ",zsh)
    print("zch: ",abs(zch),cmath.phase(zch)*(180/math.pi)," -> ",zch)
    #print("zsh: ",abs(zsh)*((13800*13800)/10000000)," -> ",zsh)
    print("Vspah:",abs(vspach),cmath.phase(vspach)*(180/math.pi)," -> ",vspach)
    print("Vcpah:",abs(vcpach),cmath.phase(vcpach)*(180/math.pi)," -> ",vcpach)
    #print("Vspah + Vcpah:",(vspach+vcpach))
    print("Vpac:",abs(vpac),cmath.phase(vpac)*(180/math.pi)," -> ",vpac)
    #print("\nVpac (V),polar:",cmath.polar(vpac*13800))
    vs_proj = projection_complex(vspach,vpac)
    vc_proj = projection_complex(vcpach,vpac)
    print("\nVs_proj:",abs(vs_proj))
    print("\nVc_proj:",abs(vc_proj))

    print("\nResp. Sitema:",(abs(vs_proj)/abs(vpac))*100)
    print("\nResp. Barra",c+1,":",(abs(vc_proj)/abs(vpac))*100)
    
    #zsh2 = (cmath.rect(6.0572,(78.9041*(math.pi/180))))
    #ipac2 = (cmath.rect(17.7580,(-146.8184*(math.pi/180))))
    #print("ISH: ", abs((((vpac)/zsh) + ipac)), cmath.phase((((vpac)/zsh) + ipac))*(180/math.pi))
    #print("meu vpac: ", abs(vpac), cmath.phase(vpac)*(180/math.pi))
    #print("seu vpac: ", abs(vpac2), cmath.phase(vpac2)*(180/math.pi))
    
#CÓDIGO

Sb = float(input("Insira a base de potência: "))
qntbars = int(input("Insira a quantidade de barras: "))
print("-="*30)
print("TIPOS DE BARRA:\n1 - SLACK (Referência) \n2 - PV \n3 - PQ")
print("-="*30)

bars = []

for i in range(qntbars):
    print("-="*30)
    print("BARRA",i+1)
    print("-="*30)
    new_bar(bars)

barsqnt = bars_qnt(bars)
print("-="*30)
print("LIGAÇÕES (-1 para parar)")
print("-="*30)

y = new_matriz(qntbars, qntbars)

connections = []
vh = []
frombar = 0

while(frombar != -1):
    frombar = int(input("Da barra: "))
    if frombar != -1:
        tobar = int(input("Para barra: "))
        new_connect(y,frombar,tobar,bars, connections)
print(connections)
#print(y)
soma_shunt(y,bars)   
#print("Connections: ",connections)
#bars[0][3] = -0.0487
#bars[2][3] = 0.1606
#bars[0][2] = 1.0329
e = pow(10,-15)
#print(create_jacob(bars,y))
#print(y)
NewtonRhapson(y,bars,barsqnt,e)
print(calc_impedancias(1,bars,5, connections))
harmonic = get_corrente(qntbars,Sb)
#print(connections)
harmonic_calc(harmonic,connections,qntbars,bars,vh)
compartilha(y,bars,connections,harmonic,vh)