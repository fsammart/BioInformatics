import os

import inquirer


from exercises.exercise1 import Exercise1
from exercises.exercise2 import Exercise2
from exercises.exercise2 import Type


class ClientService:
    # Main Menu
    EXERCISE1 = 'Exercise 1'
    EXERCISE2 = 'Exercise 2'
    EXERCISE3 = 'Exercise 3'
    EXERCISE4 = 'Exercise 4'
    PROTEIN = 'Protein'
    NUCLEIC = 'Nucleic'
    ONLINE = "Online"
    LOCAL = "Local"
    EXIT = 'Exit'
    BACK = 'Back'
    INPUT_FILES_DIRECTORY_EX1 = 'archives/exercise1/'
    DEFAULT_OUTPUT_EX1 = 'archives/prot_sequences/result.fas'
    INPUT_FILES_DIRECTORY_EX2 = 'archives/prot_sequences/'
    DEFAULT_OUTPUT_EX2 = 'archives/blast/blast.out'

    MAIN_MENU = [EXERCISE1, EXERCISE2, EXERCISE3, EXERCISE4, EXIT]
    BLAST_INPUT_FORMATS = [PROTEIN, NUCLEIC]
    PROCESSING_TYPES = [ONLINE, LOCAL]

    @staticmethod
    def get_menu_answer(menu, message='What do you want to do?'):
        questions = [inquirer.List('choice', message=message, choices=menu, )]
        return inquirer.prompt(questions)['choice']

    @staticmethod
    def get_output_file(default_output):
        questions = [inquirer.Path('output_file', message='What is the name and path of the output file?',
                                   default=default_output)]
        return inquirer.prompt(questions)['output_file']

    def start(self):
        answer = self.get_menu_answer(self.MAIN_MENU, message="Which exercise would you like to run?")
        if answer == self.EXERCISE1:
            chosen_input_file = self.INPUT_FILES_DIRECTORY_EX1 \
                                + self.choose_input_file_menu(self.INPUT_FILES_DIRECTORY_EX1)
            if chosen_input_file == self.BACK:
                self.start()
            else:
                output_file = self.get_output_file(self.DEFAULT_OUTPUT_EX1)
                Exercise1.run(chosen_input_file, output_file)

        elif answer == self.EXERCISE2:
            chosen_input_file = self.INPUT_FILES_DIRECTORY_EX2 \
                                + self.choose_input_file_menu(self.INPUT_FILES_DIRECTORY_EX2)
            if chosen_input_file == self.BACK:
                self.start()
            else:
                output_file = self.get_output_file(self.DEFAULT_OUTPUT_EX2)
                blast_format = self.get_menu_answer(menu=self.BLAST_INPUT_FORMATS,
                                                    message="Is your sequence protein or nucleotide?")
                if blast_format is self.PROTEIN:
                    blast_format = Type.PROT
                else:
                    blast_format = Type.NUC
                online = self.get_menu_answer(menu=self.PROCESSING_TYPES,
                                                    message="Where do you want to process?")
                if online is self.ONLINE:
                    online = True
                else:
                    online = False
                Exercise2.run(chosen_input_file, output_file, blast_format, online)

        elif answer == self.EXERCISE3:
            chosen_input_file = self.choose_input_file_menu(self.INPUT_FILES_DIRECTORY_EX1)
            if chosen_input_file == self.BACK:
                self.start()
            else:
                Exercise1.run(chosen_input_file, self.DEFAULT_OUTPUT_EX1)

        elif answer == self.EXERCISE4:
            chosen_input_file = self.choose_input_file_menu(self.INPUT_FILES_DIRECTORY_EX1)
            if chosen_input_file == self.BACK:
                self.start()
            else:
                Exercise1.run(chosen_input_file, self.DEFAULT_OUTPUT_EX1)

    def choose_input_file_menu(self, directory=INPUT_FILES_DIRECTORY_EX1, message = "Which input file do you choose?"):
        files = self.get_files_from_directory(directory)
        files.append(self.BACK)
        return self.get_menu_answer(files, message = message)

    @staticmethod
    def get_files_from_directory(directory):
        files = [name for name in os.listdir(directory)
                 if not os.path.isdir(os.path.join(directory, name))]
        files.sort()
        return files
