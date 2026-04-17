classdef Matrix
    methods(Static)
        function prod = matmul(m1, m1_rows, m1_cols, m2, m2_rows, m2_cols)
            assert(isequal(m1_cols, m2_rows));

            prod = zeros(m1_rows, m2_cols);
            for i = 1:m1_rows
                for j = 1:m2_cols
                    for k = 1:m1_cols
                        prod(i, j) = prod(i, j) + m1(i, k) * m2(k, j);
                    end
                end
            end
        end

        function [m_inv, pseudo] = matinv(m, m_rows, m_cols)
            assert(isequal(m_rows, m_cols));
            augmented_matrix = 
        end

        function m_mag = matmag(m, m_rows, m_cols)
        end

        function norm_m = matnorm(m, m_rows, m_cols)
        end

        function cholesky_rt = matsqrt(m, m_rows, m_cols)
        end

        function eigens = mateigens(m, m_rows, m_cols)
        end
        
        function sum = matadd(m1, m1_rows, m1_cols, m2, m2_rows, m2_cols)
        end

        function transpose = mattranspose(m, m_rows, m_cols)
        end
        
        function symmetric = matsymcheck(m, m_rows, m_cols)
        end
    end
end