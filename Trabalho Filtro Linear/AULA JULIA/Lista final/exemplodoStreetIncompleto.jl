




function linearprogram(A,b,c,tol.MaxIter,a,p,graf)
    (m,n) = size(A)

    Indx = 1:n
    Indy = n+1:n+m
    Inds = n+m+1:n+m+n
    IndComple = n+m+1:2*n+m

    #Inicialização
    x = ones(n) #x>=0
    y = ones(m)
    s = ones(n) #s>=0

    μ = p*x'*s/n
    t = 1
    z = [x;y;s]
    F = f(A,b,c,x,y,s,μ)
    fobj = -C'*x
    gap = x'*s

    ConvOpt = 0
    ConvInfeasible = 0
    ConvInfinity = 0
    ConvMaxIter = 0
    cvg = 0

    iter = 1
    while cvg <1
        Fk = [A'*y +s-c;A*x-b;Diagonal(x)*Diagonal(s)*ones(length(x),1) - ones(length(x),1)*μ]
        Jk = [zeros(n,n) A' eye(n); A zeros(m,m) zeros(m,n); Diagonal(s) zeros(n,m) Diagonal(x)]

        d = -Jk\Fk
        dx = d[Indx]
        dy = d[Indy]
        ds = d[Inds]

        #Passo implementado
        tx=t
        ts=t
        for i = 1:n
            if x[i] ++ d[Indx[i]] <0
                tx = min(tx, x[i]/abs.(d[Indx[i]]))
            end
            if s[i] + d[Inds[i]] < 0
                ts = min(ts, s[i]/abs.(d[Inds[i]]))
            end
        end
        tx = min(1,tx*a)
        ts = min(1,ts*a)

        x = x + tx*d[Indx]
        y = y + ts*d[Indy]
        s = s + ts*d[Inds]
        z = [z [x;y;s]]
        F = [F f(A,b,c,x,y,s,μ)]
        fobj = [fobj (-c'*s)]
        gap = [gap x'*s]
        iter = iter +1
        μ = p*x'*s/n



        #Olhando a convergencia
        if maximum(abs.(y))< 1e3/tol
            if gap[iter] < tol
                ConvOpt = 1
            end
        else
            ConvInfeasible =1
        end
        if maximum(abs.(x)) < 1e3/tol
            if gap[iter] < tol
                ConvOpt = 1
            end
        else
            ConvInfinity = 1
        end
        if iter >= MaxIter
            ConvMaxIter = 1
        end
        cvg = ConvOpt*1 + ConvInfeasible*2 + ConvInfinity*3 + ConvMaxIter*4
    end




    ## tem mais alguma coisa aqui


end
