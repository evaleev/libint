
#include <stdexcept>
#include <task.h>

using namespace std;
using namespace libint2;

LibraryTaskManager::LibraryTaskManager() : current_(-1)
{
}

LibraryTaskManager LibraryTaskManager::LTM_obj_;

LibraryTaskManager&
LibraryTaskManager::Instance()
{
  return LTM_obj_;
}

void
LibraryTaskManager::add(const std::string& task_label)
{
  tasks_citer end = tasks_.end();
  for(tasks_citer t=tasks_.begin(); t!=end; ++t) {
    if (t->label() == task_label)
      return;
  }
  tasks_.push_back(LibraryTask(task_label,SafePtr<TaskParameters>(new TaskParameters)));
  // make the first added task current
  if (tasks_.size() == 1)
    current_ = 0;
}

void
LibraryTaskManager::current(const std::string& task_label)
{
  tasks_citer end = tasks_.end();
  for(tasks_citer t=tasks_.begin(); t!=end; ++t) {
    if (t->label() == task_label) {
      current_ = t - tasks_.begin();
      return;
    }
  }
  throw std::runtime_error("LibraryTaskManager::current -- unknown task");
}

LibraryTask&
LibraryTaskManager::current()
{
  if (current_ >= 0)
    return tasks_.at(current_);
  else
    throw std::runtime_error("LibraryTaskManager::current -- current task has not been assigned");
}
