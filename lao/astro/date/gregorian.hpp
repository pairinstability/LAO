/// the design of this is such that the only exposed methods are `GregDay`, `GregMonth`, `GregYear`.
/// these are used to construct a gregorian date.
/// TODO: would be nice to be able to construct a month with a string like GregMonth("January") but its a pain right now

#ifndef LAO_ASTRO_DATE_GREGORIAN_H_
#define LAO_ASTRO_DATE_GREGORIAN_H_

#include <stdexcept>
#include <string>

namespace lao {
namespace astro {

    namespace detail {
        template <typename T, T Min, T Max>
        class RangeValidator {
        public:
            explicit RangeValidator(T value, const std::string& error_message)
            {
                validate(value, error_message);
                m_value = value;
            }
            using value_type = T;

            operator value_type() const
            {
                return m_value;
            }

            virtual ~RangeValidator() {};

        private:
            T m_value;

            void validate(T value, const std::string& error_message) const
            {
                if (value < Min || value > Max) {
                    throw std::out_of_range(error_message);
                }
            }
        };
    }; // namespace detail

    /// @brief represents a gregorian day.
    /// @details possible integral range is 1 to 31.
    class GregDay : public detail::RangeValidator<unsigned int, 1, 31> {
    public:
        explicit GregDay(unsigned int day)
            : RangeValidator(day, "Day must be in the range 1-31")
        {
        }
    };

    /// @brief represents a gregorian month.
    /// @details possible integral range is 1 to 21.
    class GregMonth : public detail::RangeValidator<unsigned int, 1, 12> {
    public:
        explicit GregMonth(unsigned int month)
            : RangeValidator(month, "Month must be in the range 1-12")
        {
        }
    };

    /// @brief represents a gregorian year.
    /// @details possible integral change is 1000 to 9999.
    class GregYear : public detail::RangeValidator<unsigned int, 1000, 9999> {
    public:
        explicit GregYear(unsigned int year)
            : RangeValidator(year, "Year must be in the range 1000-9999")
        {
        }
    };

} // namespace astro
} // namespace lao

#endif // LAO_ASTRO_DATE_GREGORIAN_H_