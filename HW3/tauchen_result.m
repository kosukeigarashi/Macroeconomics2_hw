N = 5; % Number of Grid Points, the number of potential realizations of z.
mu = 0; % Mean
rho = 0.9; % AR(1) Coefficient
sigma = 0.1; % Standard Deviation
m = 1; % Number of Standard Deviations


[Z,Zprob] = tauchen(N,mu,rho,sigma,m);

% 計算結果の表示
disp('状態空間のグリッド点 (Z):');
disp(Z);

disp('遷移確率行列 (P):');
disp(Zprob);

% 各行の遷移確率が正規化されていることを確認
disp('各行の遷移確率の合計 (すべて1であるべき):');
disp(sum(Zprob, 2));


% 定常分布の計算
[V, D] = eig(Zprob');   % Pの転置の固有値問題を解く
[~, idx] = max(abs(diag(D))); % 固有値1に対応する固有ベクトルを探す
pi = V(:, idx);     % 固有ベクトルを取得
pi = pi / sum(pi);  % 確率の和が1になるように正規化

% 定常分布の表示
disp('定常分布 (stationary distribution):');
disp(pi);

% 定常分布のプロット
figure;
bar(Z, pi, 'FaceColor', [0.2, 0.6, 0.8]);
xlabel('State space z_t');
ylabel('Stationary distribution pi');
title('Stationary distribution of AR(1) process by Tauchen method');
grid on;