#ifndef LAO_ASTRO_DATE_EPOCH_H_
#define LAO_ASTRO_DATE_EPOCH_H_

#include <lao/astro/date/gregorian.hpp>

namespace lao {
namespace astro {

    /// @brief Epoch class represents a moment in time.
    /// @details it uses the julian date format, supporting a julian date (JD), modified julian date (MJD)
    /// or a modified julian date 2000 (MJD).
    class Epoch {
    private:
        // JD is defined as days starting 4713BC, 1200 Jan 1.
        // MJD is defined as days starting 1858, 0000 Nov 17,
        // MJD2000 is defined as days starting 2000, 0000 Jan 1.
        enum Variant {
            JD_format,
            MJD_format,
            MJD2000_format
        };

        // use MJD2000 as the underlying date system for consistency
        double m_mjd2000;

        /// @brief converts JD to MJD2000.
        /// @details see https://en.wikipedia.org/wiki/Julian_day#Variants.
        /// @param epoch_date the epoch as a julian date.
        /// @returns the epoch as MJD2000.
        inline double JDToMJD2000(const double& epoch_date) const
        {
            return (epoch_date - 2451544.5);
        };

        /// @brief converts MJD to MJD2000.
        /// @details see https://en.wikipedia.org/wiki/Julian_day#Variants.
        /// @param epoch_date the epoch as a modified julian date.
        /// @returns the epoch as MJD2000.
        inline double MJDToMJD2000(const double& epoch_date) const
        {
            return (epoch_date - 51544);
        };

        /// @brief converts MJD2000 to JD.
        /// @details see https://en.wikipedia.org/wiki/Julian_day#Variants.
        /// @param epoch_date the epoch as a modified julian date 2000.
        /// @returns the epoch as JD.
        inline double MJD2000ToJD(const double& epoch_date) const
        {
            return (epoch_date + 2451544.5);
        };

        /// @brief converts MJD2000 to MJD.
        /// @details see https://en.wikipedia.org/wiki/Julian_day#Variants.
        /// @param epoch_date the epoch as a modified julian date 2000.
        /// @returns the epoch as MJD.
        inline double MJD2000ToMJD(const double& epoch_date) const
        {
            return (epoch_date + 51544);
        };

    public:
        /// @brief constructor for constructing an epoch from a julian date.
        /// @param epoch_date the epoch date as a julian date.
        /// @param epoch_type the type of julian date, one of Variant (JD, MJD, or MJD2000).
        Epoch(const double& epoch_date = 0, Variant epoch_type = MJD2000_format)
            : m_mjd2000(epoch_date)
        {

            switch (epoch_type) {
            case MJD2000_format:
                break;
            case JD_format:
                m_mjd2000 = JDToMJD2000(epoch_date);
                break;
            case MJD_format:
                m_mjd2000 = MJDToMJD2000(epoch_date);
                break;
            }
        };

        /// @brief constructor for constructing an epoch from a gregorian date.
        /// @details as the Epoch class uses a julian date system this constructor converts from gregorian to
        /// julian using https://en.wikipedia.org/wiki/Julian_day#Julian_day_number_calculation.
        /// @param day a gregorian day.
        /// @param month a gregorian month.
        /// @param year a gregorian year.
        Epoch(const GregDay& day, const GregMonth& month, const GregYear& year)
        {
            m_mjd2000 = (1461.0 * (year + 4800.0 + (month - 14.0) / 12.0)) / 4.0
                + (367.0 * (month - 2.0 - 12.0 * ((month - 14.0) / 12.0))) / 12.0
                - (3.0 * ((year + 4900.0 + (month - 14.0) / 12.0) / 100.0)) / 4.0
                + day - 32075.0;
        };

        /// @brief JD getter.
        /// @returns the epoch as a JD.
        double JD() const
        {
            return MJD2000ToJD(m_mjd2000);
        };

        /// @brief MJD getter.
        /// @returns the epoch as an MJD.
        double MJD() const
        {
            return MJD2000ToMJD(m_mjd2000);
        };

        /// @brief MJD2000 getter.
        /// @returns the epoch as an MJD2000.
        double MJD2000() const
        {
            return m_mjd2000;
        };
    };

} // namespace astro
} // namespace lao

#endif // LAO_ASTRO_DATE_EPOCH_H_