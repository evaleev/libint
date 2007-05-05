
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
  tasks_.push_back(LibraryTask(task_label,SafePtr<TaskParameters>(new TaskParameters),SafePtr<TaskExternSymbols>(new TaskExternSymbols)));
  // make the first added task current
  if (tasks_.size() == 1)
    current_ = 0;
}

LibraryTaskManager::TasksCIter
LibraryTaskManager::find(const std::string& task_label) const
{
  tasks_citer end = tasks_.end();
  for(tasks_citer t=tasks_.begin(); t!=end; ++t) {
    if (t->label() == task_label)
      return t;
  }
  throw ProgrammingError("LibraryTaskManager::find() -- the task not found");
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

////

void
TaskExternSymbols::add(const SymbolList& symbols)
{
  typedef SymbolList::const_iterator citer;
  citer end = symbols.end();
  for(citer s=symbols.begin(); s!=end; ++s) {
    symbols_[*s] = true;
  }
}

void
TaskExternSymbols::add(const RRList& rrlist)
{
  typedef RRList::const_iterator citer;
  citer end = rrlist.end();
  for(citer rr=rrlist.begin(); rr!=end; ++rr) {
    rrmap_[*rr] = true;
  }
}

const TaskExternSymbols::SymbolList&
TaskExternSymbols::symbols() const
{
  symbollist_.clear();
  typedef Symbols::const_iterator citer;
  citer end = symbols_.end();
  for(citer s=symbols_.begin(); s!=end; ++s) {
    symbollist_.push_back(s->first);
  }
  return symbollist_;
}

bool
TaskExternSymbols::find(const RRid& rrid) const
{
  typedef RRmap::const_iterator citer;
  citer found = rrmap_.find(rrid);
  return (found != rrmap_.end());
}

