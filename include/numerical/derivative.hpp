#ifndef DERIVATIVE_HPP_INCLUDED
#define DERIVATIVE_HPP_INCLUDED

#include <cmath>
#include <functional>
#include <optional>
#include <stdexcept>

namespace numerical_algorithms
{
    namespace derivative
    {

        // enum class that specifies the type of method used to calculate the derivative of a given function
        enum class SelectedMethod
        {
            Forward,
            Backward,
            Central
        };

        class dtobj
        {
            public:
                /********************************* CONSTRUCTORS **********************************/
                
                // default constructor
                dtobj() = default;

                // constructor with member initializer list
                dtobj(double eval_point) noexcept : _EvalPoint(eval_point) { }

                // copy constructor
                dtobj(const dtobj &) = default;

                /***************************  METHODS OF THE CLASS *******************************/
                
                /*============================== SECURITY CHECKS ================================*/

                /* Check whether the evaluation point has been initialized correctly or not */
                bool isInitializedFinite() const
                {
                    // The order of operation is the following:
                    // 1. Check if the evaluation point has been initialized
                    // 2. Check if the evaluation point contains a finite value or not
                    return _EvalPoint.has_value() && std::isfinite(_EvalPoint.value());
                }

                /*===================================== SETTER =================================*/
                void SetEvalPoint(const double eval_point)
                {
                    if(std::isfinite(eval_point))
                        _EvalPoint = eval_point;
                    else
                        throw std::invalid_argument("Evaluation point must be a finite value");
                }


                /* Calculate the first derivative with the following methods:
                *
                *  1. Forward difference
                *  2. Backward difference
                *  3. Central difference
                * 
                */
                inline double FirstDerivative(const std::function<double(double)>& function, const double spacing, SelectedMethod method_name) const
                {
                    if(!isInitializedFinite())
                        throw std::runtime_error("Evaluation point not initialized correctly.");

                    switch(method_name)
                    {
                        case SelectedMethod::Forward:
                        {
                            return ForwardDifference(function, _EvalPoint.value(), spacing);
                        }
                        case SelectedMethod::Backward:
                        {
                            return BackwardDifference(function, _EvalPoint.value(), spacing);
                        }
                        case SelectedMethod::Central:
                        {
                            return CentralDifference(function, _EvalPoint.value(), spacing);
                        }
                        default:
                        {
                            throw std::runtime_error("Method not implemented!");
                        }
                    }
                }

                inline double FourthOrderApprox(const std::function<double(double)>& function, const double spacing) const
                {
                    if(!isInitializedFinite())
                        throw std::runtime_error("Evaluation point not initialized correctly.");

                    const double f_xm2h = function(_EvalPoint.value() - 2*spacing); // f(x - 2h)
                    const double f_xmh = function(_EvalPoint.value() - spacing); // f(x - h)
                    const double f_xph = function(_EvalPoint.value() + spacing); // f(x + h)
                    const double f_xp2h = function(_EvalPoint.value() + 2*spacing); // f(x + 2h)
                    
                    if(std::isfinite(spacing) && spacing > 0 && std::isfinite(f_xm2h) && std::isfinite(f_xmh) && std::isfinite(f_xph) && std::isfinite(f_xp2h))
                        return (f_xm2h - 8*f_xmh + 8*f_xph - f_xp2h)/(12*spacing);
                    else
                        throw std::runtime_error("Divergence detected. NAN or INFINITY detected.");
                }

                inline double SecondDerivative(const std::function<double(double)>& function, const double spacing, SelectedMethod method_name) const
                {
                    if(!isInitializedFinite())
                        throw std::runtime_error("Evaluation point not initialized correctly.");

                    switch(method_name)
                    {
                        case SelectedMethod::Forward:
                        {
                            return ForwardDifferenceSecondDerivative(function, _EvalPoint.value(), spacing);
                        }
                        case SelectedMethod::Backward:
                        {
                            return BackwardDifferenceSecondDerivative(function, _EvalPoint.value(), spacing);
                        }
                        case SelectedMethod::Central:
                        {
                            return CentralDifferenceSecondDerivative(function, _EvalPoint.value(), spacing);
                        }
                        default:
                        {
                            throw std::runtime_error("Method not implemented!");
                        }
                    }
                }

            private:
                std::optional<double> _EvalPoint;

                /*  Methods to calculate the first derivative  */

                inline double ForwardDifference(const std::function<double(double)>& function, const double eval_point, const double spacing) const
                {
                    const double f_x = function(eval_point); // f(x)
                    const double f_xph = function(eval_point + spacing); // f(x + h)

                    if(std::isfinite(spacing) && spacing > 0 && std::isfinite(f_x) && std::isfinite(f_xph))
                        return (f_xph - f_x)/spacing;
                    else
                        throw std::runtime_error("Divergence detected. NAN or INFINITY detected.");
                }

                inline double BackwardDifference(const std::function<double(double)>& function, const double eval_point, const double spacing) const
                {
                    const double f_x = function(eval_point); // f(x)
                    const double f_xmh = function(eval_point - spacing); // f(x - h)

                    if(std::isfinite(spacing) && spacing > 0 && std::isfinite(f_x) && std::isfinite(f_xmh))
                        return (f_x - f_xmh)/spacing;
                    else
                        throw std::runtime_error("Divergence detected. NAN or INFINITY detected.");
                }

                inline double CentralDifference(const std::function<double(double)>& function, const double eval_point, const double spacing) const
                {
                    const double f_xph = function(eval_point + spacing); // f(x + h)
                    const double f_xmh = function(eval_point - spacing); // f(x - h)

                    if(std::isfinite(spacing) && spacing > 0 && std::isfinite(f_xph) && std::isfinite(f_xmh))
                        return (f_xph - f_xmh)/(2*spacing);
                    else
                        throw std::runtime_error("Divergence detected. NAN or INFINITY detected.");
                }

                /*  Methods to calculate the second derivative  */

                inline double ForwardDifferenceSecondDerivative(const std::function<double(double)>& function, const double eval_point, const double spacing) const
                {
                    if(!std::isfinite(spacing) || spacing < 0)
                        throw std::runtime_error("Invalid spacing provided.");
                    
                    const double f_x = function(eval_point); // f(x)
                    const double f_xph = function(eval_point + spacing); // f(x + h)
                    const double f_xp2h = function(eval_point + 2*spacing); // f(x + 2h)

                    if(std::isfinite(f_x) && std::isfinite(f_xph) && std::isfinite(f_xp2h))
                        return (f_xp2h - 2 * f_xph + f_x) / (spacing * spacing);
                    else
                        throw std::runtime_error("Divergence detected. NAN or INFINITY detected.");
                }

                inline double BackwardDifferenceSecondDerivative(const std::function<double(double)>& function, const double eval_point, const double spacing) const
                {
                    if(!std::isfinite(spacing) || spacing < 0)
                        throw std::runtime_error("Invalid spacing provided.");

                    const double f_x = function(eval_point); // f(x)
                    const double f_xmh = function(eval_point - spacing); // f(x - h)
                    const double f_xm2h = function(eval_point - 2*spacing); // f(x - 2h)

                    if(std::isfinite(f_x) && std::isfinite(f_xmh) && std::isfinite(f_xm2h))
                        return (f_x - 2 * f_xmh + f_xm2h) / (spacing * spacing);
                    else
                        throw std::runtime_error("Divergence detected. NAN or INFINITY detected.");                    
                }

                inline double CentralDifferenceSecondDerivative(const std::function<double(double)>& function, const double eval_point, const double spacing) const
                {
                    if(!std::isfinite(spacing) || spacing < 0)
                        throw std::runtime_error("Invalid spacing provided.");

                    const double f_x = function(eval_point); // f(x)
                    const double f_xmh = function(eval_point - spacing); // f(x - h)
                    const double f_xph = function(eval_point + spacing); // f(x + h)

                    if(std::isfinite(f_x) && std::isfinite(f_xph) && std::isfinite(f_xmh))
                        return (f_xph - 2 * f_x + f_xmh) / (spacing * spacing);
                    else
                        throw std::runtime_error("Divergence detected. NAN or INFINITY detected.");
                }

        }; // class dtobj

    } // namespace derivative

} // namespace numerical_algorithms


#endif // DERIVATIVE_HPP_INCLUDED
