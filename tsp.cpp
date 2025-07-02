#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <random>
#include <chrono>
using namespace std;

const double EPS = 1e-8;

class SimplexMethod {
public:
    int n;
    vector<vector<double>> costs;
    vector<double> c;
    vector<vector<double>> A;
    vector<double> b;
    vector<int> basis;
    vector<vector<double>> table;
    int num_slack_vars;
    vector<vector<int>> subtours;
    vector<int> start_basis;
    vector<vector<double>> added_gomory_cuts; // Храним добавленные ограничения
    void initialize() {
        int total_vars = n * n - n;  // Теперь переменных n² - n
        c.resize(total_vars, 0.0);

        // Заполняем вектор стоимостей c (без диагонали)
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                if (i != j) {
                    int k = i * (n - 1) + (j > i ? j - 1 : j);  // Новая нумерация
                    c[k] = costs[i][j];
                }
            }
        }

        int num_constraints = 2 * n - 1;
        A = vector<vector<double>>(num_constraints, vector<double>(total_vars, 0.0));
        b = vector<double>(num_constraints, 1.0);
        // 1. Ограничения на выход из каждого города (sum_j x_{i,j} = 1)
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                if (i == j) continue;  // Пропускаем диагональные элементы

                // Линейный индекс для x_{i,j}
                int k = i * (n - 1) + (j < i ? j : j - 1);
                A[i][k] = 1.0;
            }
        }

        // 2. Ограничения на вход в каждый город (sum_i x_{i,j} = 1)
        for (int j = 0; j < n - 1; ++j) {
            for (int i = 0; i < n; ++i) {
                if (i == j) continue;  // Пропускаем диагональные элементы

                // Тот же линейный индекс, что и выше
                int k = i * (n - 1) + (j < i ? j : j - 1);
                A[n + j][k] = 1.0;
            }
        }

        basis.clear();
        // 1. Добавляем первые (n-1) переменных с шагом (n-1)
        for (int i = 1; i <= n; ++i) {
            basis.push_back(i * (n - 1));
        }

        // 2. Добавляем следующие n переменных с шагом n (кроме уже добавленных)
        for (int i = 1; i <= n; ++i) {
            int var = (i * n) % (n * n - n);
            bool is_duplicate = false;
            for (int existing_var : basis) {
                if (existing_var == var || var==0) {
                    is_duplicate = true;
                    break;
                }
            }
            // Проверяем, не было ли это число уже добавлено на первом этапе
            if (!is_duplicate) {
                basis.push_back(var);
            }
        }

        // 3. Добавляем (n-2) переменную (если n > 2)
        if (n > 2) {
            basis.push_back(n - 2);
        }
        // Убедимся, что в базисе ровно 2n-1 переменных (как требуется в TSP)
        while (basis.size() < 2 * n - 1) {
            // Добавляем недостающие переменные по возрастанию, исключая уже имеющиеся
            for (int k = 1; k <= n * (n - 1); ++k) {
                if (find(basis.begin(), basis.end(), k) == basis.end()) {
                    basis.push_back(k);
                    if (basis.size() == 2 * n - 1) break;
                }
            }
        }
        start_basis = basis;
        
        create_table();
    }

    void create_table() {
        int total_vars = n * n - n + num_slack_vars;
        int num_rows = A.size();
        table = vector<vector<double>>(num_rows + 1, vector<double>(total_vars + 1, 0.0));

        for (int i = 0; i < num_rows; ++i) {
            for (int j = 0; j < total_vars; ++j) {
                table[i][j] = A[i][j];
            }
            table[i][total_vars] = b[i];
        }

        for (int j = 0; j < total_vars; ++j) {
            table[num_rows][j] = c[j];
        }
        table[num_rows][total_vars] = 0.0;
        //std::cout << "Current basis variables: ";
        //for (int var : basis) {
         //   std::cout << var << " ";
        //}
        //std::cout << std::endl;
        print_simplex_table();
        //приводим к каноничному виду
        for (int i = 0; i < basis.size(); ++i) {
            int var = basis[i]-1;
            double pivot_val = table[i][var];
            if (abs(pivot_val) < EPS) continue;

            for (int j = 0; j <= total_vars; ++j) {
                table[i][j] /= pivot_val;
            }

            for (int k = 0; k <= num_rows; ++k) {
                if (k == i) continue;
                double coeff = table[k][var];
                if (abs(coeff) > EPS) {
                    for (int j = 0; j <= total_vars; ++j) {
                        table[k][j] -= coeff * table[i][j];
                    }
                }
            }
        }
        
        print_simplex_table();
    }

    void pivot(int pivot_row, int pivot_col) {
        int total_vars = table[0].size() - 1;
        double pivot_val = table[pivot_row][pivot_col];
        if (abs(pivot_val) < EPS) return;

        for (int j = 0; j <= total_vars; ++j) {
            table[pivot_row][j] /= pivot_val;
        }

        for (int i = 0; i < table.size(); ++i) {
            if (i == pivot_row) continue;
            double factor = table[i][pivot_col];
            if (abs(factor) > EPS) {
                for (int j = 0; j <= total_vars; ++j) {
                    table[i][j] -= factor * table[pivot_row][j];
                }
            }
        }

        basis[pivot_row] = pivot_col+1;
    }

    vector<vector<double>> get_current_solution() {
        vector<vector<double>> sol(n, vector<double>(n, 0.0));
        int total_vars = table[0].size() - 1;

        for (int i = 0; i < basis.size(); ++i) {
            int var = basis[i]-1;
            if (var < n * (n-1)) {
                var += var / n+1;
                int r = var / n;
                int c = var % n;
                sol[r][c] = table[i][total_vars];
            }
        }
        //for (auto& row : sol) {
        //    for (double val : row) {
        //        cout << val << " ";
        //    }
         //   cout << "\n";
        //}
        return sol;
    }

    vector<vector<int>> find_subtours(vector<vector<double>>& sol) {
        vector<int> next(n, -1);
        for (int i = 0; i < n; ++i) {
            double max_val = -1.0;
            for (int j = 0; j < n; ++j) {
                if (sol[i][j] > max_val) {
                    max_val = sol[i][j];
                    next[i] = j;
                }
            }
        }

        vector<int> color(n, 0);
        vector<vector<int>> cycles;

        for (int i = 0; i < n; ++i) {
            if (color[i] == 0) {
                vector<int> path;
                int cur = i;
                while (color[cur] != 2) {
                    if (color[cur] == 1) {
                        auto it = find(path.begin(), path.end(), cur);
                        if (it != path.end()) {
                            cycles.push_back(vector<int>(it, path.end()));
                        }
                        break;
                    }
                    color[cur] = 1;
                    path.push_back(cur);
                    cur = next[cur];
                }
                for (int node : path) {
                    color[node] = 2;
                }
            }
        }

        return cycles;
    }
    int get_index(int i, int j) {
        if (i == j) return -1;  // Диагональ исключена
        return i * (n - 1) + (j < i ? j : j - 1);
    }
    int pair_to_index(int i, int j, int n) {
        if (i == j) return -1;  // Исключаем диагональ (i == j)
        return i * (n - 1) + (j < i ? j : j - 1);
    }
    void add_subtour_constraint(vector<int>& subtour) {
        //cout << "----------adding constrains------------ "<< "\n";
        int total_vars = n * (n-1) + num_slack_vars;
        vector<double> new_row(total_vars, 0.0);
        int S_size = subtour.size();

        for (int i = 0; i < S_size-1; ++i) {
            int idx = pair_to_index( subtour[i], subtour[i + 1], n);
            if (idx != -1) new_row[idx] = 1.0;
        }
        int idx = pair_to_index(subtour[S_size - 1], subtour[0], n);
        if (idx != -1) new_row[idx] = 1.0;
        for (auto& row : A) {
            row.push_back(0.0);
        }
        new_row.push_back(1.0);
        A.push_back(new_row);
        b.push_back(S_size - 1);
        c.push_back(0.0);
        num_slack_vars++;
        basis = start_basis;
        basis.push_back(total_vars+1);
        start_basis = basis ;
        create_table();
    }
    

    void add_gomory_cuts(const vector<int>& fractional_rows_indices) {
        if (fractional_rows_indices.empty()) return;

        //cout << "----------adding Gomory cuts------------\n";

        int total_vars = n * (n - 1) + num_slack_vars;
        int num_new_rows = 0;

        // Сначала подготовим все новые строки
        vector<vector<double>> new_rows;
        vector<double> new_rhs;

        for (int row_idx : fractional_rows_indices) {
            vector<double> fractional_row = table[row_idx];
            vector<double> new_row(total_vars, 0.0);

            // Вычисляем дробные части
            for (int i = 0; i < total_vars; ++i) {
                double aaa = round(fractional_row[i] * 1e10)/ 1e10;
                if (aaa == -0)
                {
                    aaa = 0;
                }
                new_row[i] =  floor(aaa)-fractional_row[i];
            }

            double rhs_fractional = fractional_row.back() - floor((fractional_row.back() * 1e10) / 1e10);

            // Проверяем, не добавляли ли уже такое ограничение
            bool is_duplicate = false;

            for (const auto& added_cut : added_gomory_cuts)
                if (equal(
                    added_cut.begin(), added_cut.end(),  // Первый диапазон
                    new_row.begin(),                     // Начало второго диапазона
                    [](double a, double b) { return std::abs(a - b) < 1e-10; })) {
                    is_duplicate = true;
                    break;
                }
            

            if (!is_duplicate) {
                new_rows.push_back(new_row);
                new_rhs.push_back(-rhs_fractional);
                added_gomory_cuts.push_back(new_row); // Запоминаем добавленный cut
                num_new_rows++;
            }
        }

        if (new_rows.empty()) {
            cout << "All cuts are duplicates - nothing to add\n";
            int rows = table.size();
            int cols = (rows > 0) ? table[0].size() : 0;
            cout << "\n" << "----------------printing b------------------" << "\n";
            for (int i = 0; i < rows; ++i) {
                if (abs(table[i][cols - 1]) < EPS) {
                    cout << "0\t";  // Нули для читаемости
                }
                else {
                    cout << table[i][cols - 1] << "\t";
                }
            }
            return;
        }

        // Расширяем матрицу A для новых slack-переменных
        for (auto& row : A) {
            row.resize(total_vars  + num_new_rows, 0.0);
        }
        size_t non_slack_count = start_basis.size() - num_slack_vars+1;
        if (non_slack_count > basis.size()) {
            non_slack_count = basis.size(); // Защита от выхода за границы
        }

        for (size_t i = 0; i < start_basis.size(); ++i) {
            basis[i] = start_basis[i];
        }
        start_basis = basis;
        // 2. Заполняем оставшуюся часть базиса последовательными номерами
        

        // Добавляем новые ограничения
        for (int k = 0; k < num_new_rows; ++k) {
            vector<double> complete_row(total_vars  + num_new_rows, 0.0);
            new_rows[k].resize(total_vars + num_new_rows, 0);
            // Копируем коэффициенты
            copy(new_rows[k].begin(), new_rows[k].end(), complete_row.begin());

            // Устанавливаем slack-переменную
            complete_row[total_vars  + k] = 1.0;

            A.push_back(complete_row);
            b.push_back(new_rhs[k]);
            c.push_back(0.0);
            basis.push_back(total_vars  + k+1);
        }
        start_basis = basis;
        num_slack_vars += num_new_rows;
        create_table();

        //cout << "Added " << num_new_rows << " new Gomory cuts\n";
        }
    void table_to_AB() {
        // Проверка размерности таблицы
        if (table.empty() || table[0].size() != n * (n - 1) + num_slack_vars + 1) {
            throw runtime_error("Некорректные размеры симплекс-таблицы");
        }

        int num_constraints = 2 * n - 1;
        int total_vars = n * (n - 1) + num_slack_vars;

        // Инициализация матрицы ограничений и правых частей
        A.assign(num_constraints, vector<double>(total_vars, 0.0));
        b.resize(num_constraints);

        // Перенос ограничений (исключая целевую строку)
        for (int i = 0; i < num_constraints; ++i) {
            if (i >= static_cast<int>(table.size()) - 1) {
                throw runtime_error("Недостаточно строк в таблице для всех ограничений");
            }

            // Копирование коэффициентов
            copy(table[i].begin(), table[i].begin() + total_vars, A[i].begin());

            // Копирование правой части
            b[i] = table[i].back();
        }

        // Перенос целевой функции в вектор c
        int last_row = table.size() - 1;
        c.assign(total_vars, 0.0);
        copy(table[last_row].begin(), table[last_row].begin() + total_vars, c.begin());



    };
    pair<int, int> get_ij_from_index(int k) {
        if (k < 0 || k >= n * (n - 1)) {
            throw out_of_range("Invalid variable index");
        }

        int i = k / n;  // Определяем строку
        int j = k % n;  // Позиция в строке

        if (i == j) {
            throw std::out_of_range("Индекс соответствует диагональному элементу");
        }

        return make_pair(i, j);
    }
    void print_simplex_table() {
        int rows = table.size();
        int cols = (rows > 0) ? table[0].size() : 0;
        /*cout << "\n" << "----------------printing b------------------" << "\n";
        for (int i = 0; i < rows; ++i) {
               if (abs(table[i][cols-1]) < EPS) {
                    cout << "0\t";  // Нули для читаемости
               }
               else {
                    cout << table[i][cols - 1] << "\t";
               }
        }*/
        if (false) {
            // Определяем размеры таблицы
            

            // Вывод заголовков столбцов (переменные)
            cout << "        ";
            for (int j = 0; j < cols - 1; ++j) {
                if (j < n * (n - 1)) {
                    // Основные переменные x_{i,j}
                    cout << "x_" << j << "\t";
                }
                else {
                    // Slack-переменные
                    cout << "s_" << (j - n * (n - 1)) << "\t";
                }
            }
            cout << "RHS" << endl;  // Столбец правой части

            // Вывод строк таблицы
            for (int i = 0; i < rows; ++i) {
                // Вывод базисной переменной для строки
                if (i < basis.size()) {
                    int var = basis[i] - 1;
                    if (var < n * (n - 1)) {

                        cout << "x_" << var << ":\t";
                    }
                    else {
                        cout << "s_" << (var - n * (n - 1)) << ":\t";
                    }
                }
                else {
                    cout << "   \t";  // Для строки целевой функции
                }

                // Вывод коэффициентов
                for (int j = 0; j < cols; ++j) {
                    if (abs(table[i][j]) < EPS) {
                        cout << "0\t";  // Нули для читаемости
                    }
                    else {
                        cout << table[i][j] << "\t";
                    }
                }
                cout << endl;

                // Разделитель после последней строки ограничений
                if (i == rows - 2) {
                    cout << string(10 + 8 * cols, '-') << endl;
                }
            }
            cout << endl;
        }
        
    }
    bool handle_negative_rhs() {
        bool has_negative = false;
        const int rhs_col = table[0].size() - 1;
        const int last_row = table.size() - 1;

        // 1. Проверяем наличие отрицательных RHS (кроме последней строки)
        for (int i = 0; i < last_row; ++i) {
            if (table[i][rhs_col] < -EPS) {
                has_negative = true;
                break;
            }
        }

        if (!has_negative) {
            return false; // Отрицательных элементов нет
        }

        // 2. Обработка отрицательных RHS
        while (true) {
            int pivot_row = -1;
            int pivot_col = -1;
            double min_ratio = numeric_limits<double>::max();

            // Ищем строку с отрицательным RHS
            for (int i = 0; i < last_row; ++i) {
                if (table[i][rhs_col] < -EPS) {
                    // Ищем подходящий столбец для pivot
                    for (int j = 0; j < rhs_col; ++j) {
                        // Ищем отрицательный коэффициент в небазисной переменной
                        if (table[i][j] < -EPS && !is_basic_variable(j)) {
                            double ratio = abs(table[i][rhs_col] / table[i][j]);
                            if (ratio < min_ratio) {
                                min_ratio = ratio;
                                pivot_row = i;
                                pivot_col = j;
                            }
                        }
                    }
                }
            }

            if (pivot_row == -1) {
                // Не найдено допустимого pivot - задача не имеет допустимого решения
                cout<<"No solution"<<"\n";
                break;
            }

            // Выполняем pivot операцию
            pivot(pivot_row, pivot_col);

            // Проверяем, остались ли отрицательные RHS
            has_negative = false;
            for (int i = 0; i < last_row; ++i) {
                if (table[i][rhs_col] < -EPS) {
                    has_negative = true;
                    break;
                }
            }

            if (!has_negative) {
                break; // Все RHS стали неотрицательными
            }
        }
        
        return true;
    }
    vector<int> find_fractional_edges() {
        vector<int> fractional_edges;
        const int rhs_col = table[0].size() - 1;
        const int edge_vars = n * (n - 1); // Количество переменных-ребер

        for (size_t i = 0; i < basis.size(); ++i) {
            int var = basis[i];
            if (var < edge_vars) { // Проверяем только ребра
                double value = table[i][rhs_col];
                if (abs(round(value) - value) > EPS) {
                    //fractional_edges.push_back(var);
                    fractional_edges.push_back(i);
                }
            }
        }

        return fractional_edges;
    }
    bool is_basic_variable(int col) const {
        return find(basis.begin(), basis.end(), col+1) != basis.end();
    }
    SimplexMethod(vector<vector<double>> _costs) : costs(_costs) {
        n = costs.size();
        num_slack_vars = 0;
        initialize();
    }

    pair<vector<vector<double>>, double> solve(int max_iter = 1000) {
        for (int iter = 0; iter < max_iter; ++iter) {
            if (!handle_negative_rhs()) {
                int num_rows = table.size() - 1;
                int num_cols = table[0].size() - 1;
                double min_coeff = 1;
                int pivot_col = -1;

                for (int j = 0; j < num_cols; ++j) {
                    // Проверяем только коэффициенты целевой функции
                    // (все столбцы уже соответствуют допустимым переменным x_{i,j}, i!=j)
                    if (table[num_rows][j] < min_coeff) {
                        min_coeff = table[num_rows][j];
                        pivot_col = j;
                    }
                }

                if (min_coeff >= 0) break;

                double min_ratio = 10000;
                int pivot_row = -1;
                for (int i = 0; i < num_rows; ++i) {
                    if (table[i][pivot_col] > 0 and table[i][num_cols]>=0) {
                        double ratio = table[i][num_cols] / table[i][pivot_col];
                        if (ratio < min_ratio ) {
                            min_ratio = ratio;
                            pivot_row = i;
                        }
                    }
                }
                if (pivot_row == -1) {
                    break;
                }
                if (min_ratio >= 0) {
                
                    pivot(pivot_row, pivot_col);
                    
                }
                print_simplex_table();
                

                
            }
        }
            

        vector<vector<double>> sol = get_current_solution();
        vector<vector<int>> found_subtours = find_subtours(sol);
        //cout << "--------solution---------" << "\n";
        print_simplex_table();
        vector<int> fractional_edges = find_fractional_edges();
        if (fractional_edges.size() != 0) {
            add_gomory_cuts(fractional_edges);
        }
        if (fractional_edges.size() != 0) {
            return solve(max_iter - 1);
            sol = get_current_solution();
            found_subtours = find_subtours(sol);
        }
        bool has_subtours = false;
        for (auto& subtour : found_subtours) {
            if (subtour.size() < n) {
                has_subtours = true;
                //cout << "-----------has_subtours---------" << "\n";
                add_subtour_constraint(subtour);
                break;
            }
        }
        
        if (has_subtours) {
            return solve(max_iter - 1);
            sol = get_current_solution();
            found_subtours = find_subtours(sol);
        }
        
        
        double obj_value = -table.back().back();
        return make_pair(sol, obj_value);
    }
};
void RandomfillMatrix(std::vector<std::vector<double>>& matrix) {
    int N = matrix.size();
    srand(time(0));
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            matrix[i][j] = rand() % (50 - 1 + 1) + 2;  // Верхний треугольник                  // Нижний треугольник (симметрия)
        }
        matrix[i][i] = 0;
    }
}
void PrintMatrix(std::vector<std::vector<double>>& matrix) {
    for (const auto& row : matrix) {
        for (int val : row) {
            std::cout << val << " ";
        }
        std::cout << std::endl;
    }
    std::cout << "\n" << "___________________________" << "\n";
}
int main() {
    cout << numeric_limits<double>::epsilon << endl;;
    vector<vector<double>> costs = {
        { 0, 14, 16, 41, 37, 39},
        {21,  0, 42, 33, 19, 20},
        { 8, 20,  0, 27, 46, 49},
        {13, 33, 20,  0, 23,  6},
        {12, 32, 43, 24,  0,  4},
        {19, 49, 48, 47, 45,  0}
    };
    //vector<vector<double>> costs = std::vector<std::vector<double>>(10, std::vector<double>(10, 0.0));
    //RandomfillMatrix(costs);
    //PrintMatrix(costs);
    SimplexMethod solver(costs);
    std::cout << std::endl << "start timer" << std::endl;
    auto start = std::chrono::steady_clock::now();
    auto result = solver.solve(costs.size()* costs.size()*10);
    auto end = std::chrono::steady_clock::now();
    std::cout << std::endl << "end timer" << std::endl;
    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
    std::cout << std::endl << duration.count() << std::endl;
    vector<vector<double>> solution = result.first;
    double objective = result.second;

    cout << "Optimal solution:\n";
   // for (auto& row : solution) {
    //    for (double val : row) {
    //        cout << val << " ";
    //    }
    //    cout << "\n";
    //}
    cout << "\n" << "----------------" << "\n";
    vector<int> first_subtour = solver.find_subtours(solution)[0];
    for (int node : first_subtour) {
        cout << node << " ";
    }
    cout << endl;
    cout << "Total cost: " << objective << "\n";

    return 0;
}