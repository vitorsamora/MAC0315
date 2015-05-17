% Função principal e pré-processamento
function ret = simplex(A, b, c, m, n, x)
   index = FindIndex(m, n, x);
   compl = FindCompl(m, n, x);

   B = BasisMatrix(A, index, m);
   Binv = inv(B);
   ret = Simplex(A, c, m, n, x, index, compl, Binv);

   % Exibe as informações finais
   if (ret(1) == 0)
      printf("\nSolução ótima encontrada com custo: %5f\n", c' * ret(2:n+1));
   elseif(ret(1) == -1)
      printf("\nProblema ilimitado. Direção de custo -inf:\n");
   endif
   for (i = 2:n+1)
      printf("%d   %5f\n", i-1, ret(i));
   endfor
   printf("\n");
endfunction

% Método simplex revisado
function ret = Simplex(A, c, m, n, x, index, compl, Binv )
   iteration = 0;
   while (true)
      printf("\nIterando %d\n", iteration);
      indexPrint = sort(index);
      for (i = 1:m)
         printf("%d   %5f\n", indexPrint(i), x(indexPrint(i)));
      endfor

      printf("\nValor da função objetivo: %5f\n", c' * x);

      c_b = BasisVector(c, index, m);
      p = c_b' * Binv;
      c_ = c' - p * A;

      printf("\nCustos reduzidos: \n");
      for (i = 1:n-m)
         printf("%d   %5f\n", compl(i), c_(compl(i)));
      endfor

      if (c_ >= 0)
         ret = [0; x]; % solução ótima
         break;
      endif

      for (j = 1:n)
         if (c_(j) < 0)
            break;
         endif
      endfor

      u = Binv * A(:, j);

      d = NewDirec(u, index, j, m, n);

      if (u <= 0)
         ret = [-1; d]; % custo ótimo -inf
         break;
      endif

      theta = Inf;
      for (i = 1:m)
         if (u(i) > 0)
            div = x(index(i)) ./ u(i);
            if (div < theta)
               theta = div;
               imin = i;
            endif
         endif
      endfor

      printf("\nEntra na base: %d\n", j);
      printf("\nDireção:\n");
      for (i = 1:m)
         printf("%d   %5f\n", index(i), d(index(i)));
      endfor
      printf("\nTheta*: %5f\n", theta);
      printf("\nSai da base: %d\n", index(imin));

      for (i = 1:m)
         x(index(i)) = x(index(i)) - u(i) * theta;
      endfor
      x(j) = theta;

      Binv = NewBinv(Binv, u, imin, m);
      index(imin) = j;

      iteration = iteration + 1;
   endwhile
endfunction

%%%%%%%%%%%%%%% Funções auxiliares %%%%%%%%%%%%%%%

% Constrói o vetor contendo os índices das variáveis básicas
function ret = FindIndex(m, n, x)
   index = zeros(m,1);
   k = 1;
   for (i = 1:n)
      if (x(i) != 0)
         index(k) = i;
         k = k + 1;
      endif
   endfor
   ret = index;
endfunction

% Constrói o vetor contendo os índices das variáveis não básicas.
% Este vetor será usado apenas para exibir a direção das variáveis não básicas.
function ret = FindCompl(m, n, x)
   compl = zeros(n-m,1);
   k = 1;
   for (i = 1:n)
      if (x(i) == 0)
         compl(k) = i;
         k = k + 1;
      endif
   endfor
   ret = compl;
endfunction

% Constrói a matriz básica
function ret = BasisMatrix(A, index, m)
   M = zeros(m, m);
   for (i = 1:m)
      M(:, i) = A(:, index(i));
   endfor
   ret = M;
endfunction

% Constrói um vetor básico (com índices das variáveis básicas)
function ret = BasisVector(v, index, m)
  w = zeros(m, 1);
  for (i = 1:m)
     w(i) = v(index(i));
  endfor
  ret = w;
endfunction

% Constrói a matriz nversa do novo B, dada a matriz inversa do B anterior
function ret = NewBinv(Binv, u, imin, m)
   Binv(imin, :) = Binv(imin, :) / u(imin);
   for (i = 1:m)
      if (i != imin)
         Binv(i, :) = Binv(i, :) - Binv(imin, :) * u(i);
      endif
   endfor
   ret = Binv;
endfunction

% Constrói o vetor de direções
function ret = NewDirec(u, index, j, m, n)
   d = zeros(n,1);
   for (i = 1:m)
      d(index(i)) = - u(i);
   endfor
   d(j) = 1;
   ret = d;
endfunction
