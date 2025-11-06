#ifndef ROOT_FINDER_HPP_INCLUDED
#define ROOT_FINDER_HPP_INCLUDED

#include <cmath>
#include <stdexcept>

namespace numerical_algorithms
{

    namespace root_finder
    {

        /*******************************  ROOT FINDER OBJECT *************************************/
        class rfobj
        {
            public:
                /********************************* CONSTRUCTORS **********************************/
                
                // default constructor
                rfobj() = default;

                // constructor with member initializer list
                rfobj(double lower_bnd, double upper_bnd, double tol) noexcept :
                _LowerBound(lower_bnd), _UpperBound(upper_bnd), _Epsilon(tol) { }

                // copy constructor
                rfobj(const rfobj &) = default;

                /***************************  METHODS OF THE CLASS *******************************/
                
                inline double BisectionMethod(double (*function)(double)) const
                {
                    if(!check_bounds())
                        throw std::out_of_range("Invalid interval: lower bound must be less than upper bound\n");
                    if(!check_sign(function))
                        throw std::invalid_argument("The function must change its sign at the endpoints\n");

                    double left_bound = _LowerBound;
                    double right_bound = _UpperBound;

                    double midpoint = 0.0;

                    while(fabs(left_bound - right_bound) > _Epsilon)
                    {
                        midpoint = (left_bound + right_bound)*0.5;

                        if(function(left_bound)*function(midpoint) > 0)
                            left_bound = midpoint;
                        else
                            right_bound = midpoint;
                    }
                    
                    return (left_bound + right_bound)*0.5;
                }

                inline double FalsePositionMethod(double (*function)(double)) const
                {
                    if(!check_bounds())
                        throw std::out_of_range("Invalid interval: lower bound must be less than upper bound\n");
                    if(!check_sign(function))
                        throw std::invalid_argument("The function must change its sign at the endpoints\n");

                    double left_bound = _LowerBound;
                    double right_bound = _UpperBound;

                    double root_estimate = left_bound; // Initialize the root with the left bound of the interval

                    while(fabs(left_bound - right_bound) > _Epsilon)
                    {
                        root_estimate = (left_bound*function(right_bound) - right_bound*function(left_bound))/(function(right_bound) - function(left_bound));

                        if(function(left_bound)*function(root_estimate) < 0)
                            right_bound = root_estimate;
                        else
                            left_bound = root_estimate;
                    }

                    return root_estimate;
                }

                inline double SecantMethod(double (*function)(double), const unsigned int max_iter) const
                {
                    if(!check_bounds())
                        throw std::out_of_range("Invalid interval: lower bound must be less than upper bound\n");
                    if(!check_sign(function))
                        throw std::invalid_argument("The function must change its sign at the endpoints\n");
                    
                    double left_bound = _LowerBound;
                    double right_bound = _UpperBound;
                    double f_0 = function(left_bound);
                    double f_1 = function(right_bound);

                    if(f_0 == f_1)
                        throw std::invalid_argument("Attention! Division by zero detected!\n");
                    
                    unsigned int cnt_iter = 0;

                    double root_estimate = right_bound; // current estimate for the root

                    while((fabs(left_bound - right_bound) > _Epsilon))
                    {                    
                        if(cnt_iter > max_iter)
                            throw std::out_of_range("The number of iterations has exceed the maximum value\n");

                        // compute the estimate for the root
                        root_estimate = right_bound - f_1*(right_bound - left_bound)/(f_1 - f_0);
                        
                        // update everything for the next iteration
                        left_bound = right_bound;
                        f_0 = f_1;
                        right_bound = root_estimate;
                        f_1 = function(right_bound);
                        
                        // increment the number of iterations
                        cnt_iter += 1;
                    }

                    return root_estimate;
                }

                inline double NewtonRaphsonMethod(double (*function)(double), double (*derivative)(double),
                                                    double init_guess, const unsigned int max_iter) const
                {
                    if(!check_bounds())
                        throw std::out_of_range("Invalid interval: lower bound must be less than upper bound\n");
                    if(!check_sign(function))
                        throw std::invalid_argument("The function must change its sign at the endpoints\n");

                    double f_0 = function(init_guess);
                    double d_0 = derivative(init_guess);
                    double root_estimate = 0.0;

                    unsigned int cnt_iter = 0;

                    while(fabs(f_0) > _Epsilon)
                    {                    
                        if(fabs(d_0) < 1e-15)
                            throw std::runtime_error("Attention! Division by zero detected!\n");
                        
                        if(cnt_iter > max_iter)
                            throw std::out_of_range("The number of iterations has exceed the maximum value\n");
                        
                        // compute the estimate for the root
                        root_estimate = init_guess - f_0/d_0;

                        if(!std::isfinite(root_estimate))
                            throw std::runtime_error("The method diverged! NAN or INFINITY detected as value of root_estimate\n");
                        
                        if (std::fabs(root_estimate - init_guess) < _Epsilon)
                            return root_estimate; // the method has converged

                        // update everything for the next iteration
                        init_guess = root_estimate;
                        f_0 = function(init_guess);
                        d_0 = derivative(init_guess);

                        // increment the number of iterations
                        cnt_iter += 1;
                    }
                    
                    return root_estimate;
                }

            private:

                /*  The interval must be given in the mathematical form: [_LowerBound, _UpperBound]  */
                double _LowerBound;
                double _UpperBound;
                double _Epsilon; // tolerance

                /*  Private methods  */
                constexpr bool check_sign(double (*function)(double)) const noexcept
                {
                    return function(_LowerBound)*function(_UpperBound) <= 0;
                }

                constexpr bool check_bounds() const noexcept
                {
                    return _LowerBound < _UpperBound;
                }
        };

    } // namespace root_finder

} // namespace numerical_algorithms

#endif // ROOT_FINDER_HPP_INCLUDED
