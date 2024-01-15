#ifndef LOGGER_H
#define LOGGER_H

#include <sstream>
#include <iostream>
#include <chrono>
#include <ctime>

enum class TLogLevel
{
    ERROR,
    WARNING,
    DEBUG,
    INFO,
    DEBUG_TEMP
};

class LogCustom
{
public:
    class LogStream
    {
    public:
        LogStream(TLogLevel level, TLogLevel &reportingLevel, std::ostringstream &os)
            : messageLevel(level), rLevel(reportingLevel), logStream(os)
        {
        }

        ~LogStream()
        {
            if (messageLevel >= rLevel)
            {
                logStream << std::endl;
                std::cerr << logStream.str();
            }
        }

        template <typename T>
        LogStream &operator<<(const T &msg)
        {
            logStream << msg;
            return *this;
        }

    private:
        TLogLevel messageLevel;
        TLogLevel &rLevel;
        std::ostringstream &logStream;
    };

    LogCustom(TLogLevel level) : messageLevel(level) {}

    LogStream Get(TLogLevel level)
    {
        os.str("");
        os.clear();
        os << GetCurrentTime() << " " << ToString(level) << ": ";
        messageLevel = level;
        return LogStream(level, ReportingLevel(), os);
    }

    static TLogLevel &ReportingLevel()
    {
        static TLogLevel reportingLevel = TLogLevel::INFO;
        return reportingLevel;
    }

private:
    std::ostringstream os;
    TLogLevel messageLevel;

    std::string GetCurrentTime()
    {
        auto now = std::chrono::system_clock::now();
        auto now_c = std::chrono::system_clock::to_time_t(now);
        std::tm *tm = std::localtime(&now_c);

        std::ostringstream oss;
        oss << std::put_time(tm, "%d/%m/%Y %H:%M:%S");
        return oss.str();
    }

    std::string ToString(TLogLevel level)
    {
        switch (level)
        {
        case TLogLevel::ERROR:
            return "ERROR";
        case TLogLevel::WARNING:
            return "WARNING";
        case TLogLevel::INFO:
            return "INFO";
        case TLogLevel::DEBUG:
            return "DEBUG";
        default:
            return "";
        }
    }
};

#endif // LOGGER_H