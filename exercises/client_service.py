import os

import inquirer
from exercises.exercise1 import Exercise1


class ClientService:

    # Main Menu
    EXERCISE1 = 'Exercise 1'
    EXERCISE2 = 'Exercise 2'
    EXERCISE3 = 'Exercise 3'
    EXERCISE4 = 'Exercise 4'
    EXIT = 'Exit'
    BACK = 'Back'
    INPUT_FILES_DIRECTORY_EX1 = 'archives/exercise1/'
    MAIN_MENU = [EXERCISE1, EXERCISE2, EXERCISE3, EXERCISE4, EXIT]

    @staticmethod
    def get_menu_answer(menu):
        questions = [inquirer.List('choice', message="What do you want to do?", choices=menu, )]
        return inquirer.prompt(questions)['choice']

    @staticmethod
    def get_output_file():
        questions = [inquirer.Path('output_file', message='What is the name and path of the output file?', default='output/result.fas')]
        return inquirer.prompt(questions)['output_file']

    def start(self):
        answer = self.get_menu_answer(self.MAIN_MENU)
        if answer == self.EXERCISE1:
            chosen_input_file = self.INPUT_FILES_DIRECTORY_EX1 + self.choose_input_file_menu(self.INPUT_FILES_DIRECTORY_EX1)
            if chosen_input_file == self.BACK:
                self.start()
            else:
                output_file = self.get_output_file()
                Exercise1.run(chosen_input_file, output_file)

        elif answer == self.EXERCISE2:
            chosen_input_file = self.choose_input_file_menu(self.INPUT_FILES_DIRECTORY_EX1)
            if chosen_input_file == self.BACK:
                self.start()
            else:
                Exercise1.run(chosen_input_file, 'output/result1.fas')

        elif answer == self.EXERCISE3:
            chosen_input_file = self.choose_input_file_menu(self.INPUT_FILES_DIRECTORY_EX1)
            if chosen_input_file == self.BACK:
                self.start()
            else:
                Exercise1.run(chosen_input_file, 'output/result1.fas')

        elif answer == self.EXERCISE4:
            chosen_input_file = self.choose_input_file_menu(self.INPUT_FILES_DIRECTORY_EX1)
            if chosen_input_file == self.BACK:
                self.start()
            else:
                Exercise1.run(chosen_input_file, 'output/result1.fas')

    def choose_input_file_menu(self, directory=INPUT_FILES_DIRECTORY_EX1):
        files = self.get_files_from_directory(directory)
        files.append(self.BACK)
        return self.get_menu_answer(files)

    @staticmethod
    def get_files_from_directory(directory):
        files = [name for name in os.listdir(directory)
                 if not os.path.isdir(os.path.join(directory, name))]
        files.sort()
        return files
