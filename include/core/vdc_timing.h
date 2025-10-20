#ifndef VDC_TIMING_H
#define VDC_TIMING_H

#include <string>
#include <chrono>
#include <map>
#include <vector>
#include <memory>

//! @brief Timer node for hierarchical timing measurement
struct TimerNode {
    std::string name;
    std::chrono::high_resolution_clock::time_point start_time;
    double elapsed_seconds;
    bool is_running;
    TimerNode* parent;
    std::vector<std::unique_ptr<TimerNode>> children;

    TimerNode(const std::string& name, TimerNode* parent = nullptr);
    void start();
    void stop();
    double getElapsed() const;
};

//! @brief Singleton class for managing hierarchical timing measurements
class TimingStats {
public:
    static TimingStats& getInstance();

    // Prevent copying and assignment
    TimingStats(const TimingStats&) = delete;
    TimingStats& operator=(const TimingStats&) = delete;

    //! Start a timer with the given name. If parent is specified, this timer
    //! becomes a child of that parent. If parent is empty, it's a top-level timer.
    void startTimer(const std::string& name, const std::string& parent = "");

    //! Stop the timer with the given name
    void stopTimer(const std::string& name);

    //! Print the hierarchical timing report
    void printReport() const;

    //! Reset all timers (useful for testing)
    void reset();

private:
    TimingStats();
    ~TimingStats() = default;

    TimerNode* findTimer(const std::string& name);
    TimerNode* findTimerInSubtree(TimerNode* node, const std::string& name);
    void printNode(const TimerNode* node, int indent, bool is_last_child, const std::vector<bool>& ancestor_continues) const;
    std::string formatTime(double seconds) const;

    std::unique_ptr<TimerNode> root_;
    std::map<std::string, TimerNode*> timer_map_;  // For fast lookup
};

//! @brief RAII helper for automatic timer start/stop
class ScopedTimer {
public:
    ScopedTimer(const std::string& name, const std::string& parent = "");
    ~ScopedTimer();

private:
    std::string name_;
};

#endif // VDC_TIMING_H
