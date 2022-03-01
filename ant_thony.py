#!/usr/bin/env python3

"""A tool for parallel execution of tasks"""

import os
import argparse
import logging
from multiprocessing import Process, cpu_count

logging.basicConfig(
    format="[\U0001F41C-Thony] %(levelname)s: %(message)s", level=logging.DEBUG
)


class Task:
    """A task class"""

    def __init__(self, command, path="."):
        self.command = command
        self.path = path

    def run(self):
        """Runs a command in the given path"""
        os.chdir(self.path)
        os.system(self.command)


class Ant(Process):
    """Ant-Thony's buddies"""

    created = 0

    def __init__(self, tasks):
        super().__init__(name="\U0001F41C-{}".format(Ant.created+1))
        self.tasks = tasks
        logging.info("{} ready with {} tasks".format(self.name, len(self.tasks)))
        Ant.created += 1

    def run(self):
        """Runs all the assigned tasks"""
        for task in self.tasks:
            task.run()
        logging.info("{} going back to the nest".format(self.name))


class Ant_Thony:
    """Our buddy Ant-Thony"""

    def __init__(self, tasks, num_cpus=0):
        try:
            self.num_processes = int(num_cpus)
            if self.num_processes < 1:
                raise ValueError()
        except (ValueError, TypeError):
            logging.warning(
                "Number of cores has not been specified or it is incorrect. Using all available cores."
            )
            self.num_processes = cpu_count()

        logging.info("\U0001F41C-Thony will use {} cores".format(self.num_processes))

        self.tasks = tasks
        self.num_tasks = len(tasks)
        self.workers = []
        workers_tasks = [
            tasks[i :: self.num_processes] for i in range(self.num_processes)
        ]

        for i in range(self.num_processes):
            worker = Ant(workers_tasks[i])
            self.workers.append(worker)

    def release(self):
        logging.info("Swarming!")
        for ant in self.workers:
            ant.start()

        for ant in self.workers:
            ant.join()

        logging.info("{} tasks done".format(self.num_tasks))

    def go_home(self):
        for ant in self.workers:
            ant.terminate()
        logging.info("All \U0001F41C back to the nest")


if __name__ == "__main__":

    parser = argparse.ArgumentParser(prog="ant_thony")
    parser.add_argument(
        "tasks_file_name",
        help="A file containing a task for each line",
        metavar="tasks_file_name",
    )
    parser.add_argument(
        "--cores",
        "-cores",
        "-c",
        help="CPU cores to use",
        dest="cores",
        type=int,
        default=0,
    )

    args = parser.parse_args()

    with open(args.tasks_file_name) as handle:
        all_tasks = []
        for line in handle:
            if line and not line.startswith("#"):
                all_tasks.append(Task(line.rstrip(os.linesep)))

        anthony = Ant_Thony(all_tasks, args.cores)
        anthony.release()
        anthony.go_home()
