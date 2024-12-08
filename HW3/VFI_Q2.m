N = 5; % Number of Grid Points, the number of potential realizations of z.
mu = 0; % Mean
rho = 0.9; % AR(1) Coefficient
sigma = 0.1; % Standard Deviation
m = 3; % Number of Standard Deviations

[Z,Zprob] = tauchen(N,mu,rho,sigma,m);

% パラメータ設定
beta = 0.95;                % 割引因子
tau = 0.05;                 % 所得税率
alpha = 0.4;                % 資本分配率
u = @(c) log(c);            % 効用関数: log(c)
max_iter = 1000;            % 最大反復回数
tol = 1e-6;                 % 収束許容誤差

% 状態と消費の関係を表すグリッドの生成
[z_mat, z_prime_mat] = meshgrid(Z, Z); % グリッド行列化
c = z_mat - z_prime_mat;                        % 消費: 現在のz - 次期のz'

% 消費が正でない場合は効用を無効化
c(c <= 0) = NaN;

% 効用行列の計算
U = u(c); % 効用行列

% 初期化
V = zeros(N, 1);       % 初期価値関数 (ゼロベクトル)
policy = zeros(N, 1);  % 政策関数 (次期の選択)

% 価値関数反復
for iter = 1:max_iter
    V_new = zeros(N, 1);  % 更新後の価値関数
    for i = 1:N
        % 次期の価値の期待値
        EV = P(i, :) * V;

        % 現在の効用 + 割引期待値
        total_value = U(i, :) + beta * EV;

        % 最大値とそのインデックスを取得
        [V_new(i), policy(i)] = max(total_value);
    end

    % 収束判定
    if max(abs(V_new - V)) < tol
        disp(['収束しました (反復回数: ', num2str(iter), ')']);
        break;
    end
    V = V_new; % 更新
end

% 結果の表示
disp('最適価値関数:');
disp(V);

disp('最適政策関数 (次期の状態インデックス):');
disp(policy);

% 最適政策のプロット
optimal_policy = grid_z(policy); % 次期の最適状態
figure;
plot(grid_z, optimal_policy, 'LineWidth', 2);
xlabel('現在の状態 z');
ylabel('最適次期状態 z''');
title('最適政策関数');
grid on;

% 価値関数のプロット
figure;
plot(grid_z, V, 'LineWidth', 2);
xlabel('状態 z');
ylabel('価値関数 V(z)');
title('価値関数');
grid on;
