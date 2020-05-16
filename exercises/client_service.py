import os
import inquirer
from exercises.exercise1 import Exercise1
from exercises.exercise2 import Exercise2
from exercises.exercise3 import Exercise3
from exercises.exercise4 import Exercise4


class ClientService:
    # Main Menu
    EXERCISE1 = 'Exercise 1'
    EXERCISE2 = 'Exercise 2'
    EXERCISE3 = 'Exercise 3'
    EXERCISE4 = 'Exercise 4'
    EXERCISE5 = 'Exercise 5'
    EXERCISE6 = 'Exercise 6'
    PROTEIN = 'Protein'
    NUCLEIC = 'Nucleic'
    ONLINE = "Online"
    LOCAL = "Local"
    EXIT = 'Exit'
    BACK = 'Back'

    INPUT_FILES_DIRECTORY_EX1 = 'archives/gb_files/'
    DEFAULT_OUTPUT_EX1 = 'archives/prot_sequences/result.fas'
    INPUT_FILES_DIRECTORY_EX2 = 'archives/prot_sequences/'
    INPUT_FILES_DIRECTORY_EX3 = 'archives/msa/'
    INPUT_FILES_DIRECTORY_EX4 = 'archives/blast/'
    DEFAULT_OUTPUT_EX2 = 'archives/blast/blast.out'
    DEFAULT_OUTPUT_EX3 = 'archives/clustal/msa.out'
    DEFAULT_OUTPUT_EX4 = 'archives/blast/blast_pattern_hits'
    CMD_BLAST_LINUX = "/usr/bin/blastp"
    CMD_BLAST_MACOS = "/usr/local/ncbi/blast/bin/blastp"

    EX4_RESULT_SIZE = 10

    STARTED_PROCESSING = 'Processing Started'
    ENDING_PROCESSING = 'Processing Ended'

    MAIN_MENU = [EXERCISE1, EXERCISE2, EXERCISE3,EXERCISE4, EXIT]
    PROCESSING_TYPES = [ONLINE, LOCAL]

    @staticmethod
    def get_menu_answer(menu, message='What do you want to do?'):
        questions = [inquirer.List('choice', message=message, choices=menu, )]
        return inquirer.prompt(questions)['choice']

    @staticmethod
    def get_output_file(default_output, message='What is the name and path of the output file?'):
        questions = [inquirer.Path('output_file', message=message,
                                   default=default_output)]
        return inquirer.prompt(questions)['output_file']

    def start(self):
        answer = self.get_menu_answer(self.MAIN_MENU, message="Which exercise would you like to run?")

        if answer == self.EXERCISE1:
            chosen_input_file = self.choose_input_file_menu(self.INPUT_FILES_DIRECTORY_EX1)

            if chosen_input_file == self.BACK:
                self.start()
            else:
                chosen_input_file = self.INPUT_FILES_DIRECTORY_EX1 + chosen_input_file
                output_file = self.get_output_file(self.DEFAULT_OUTPUT_EX1)
                print(self.STARTED_PROCESSING)
                Exercise1.run(chosen_input_file, output_file)
                print(self.ENDING_PROCESSING)

        elif answer == self.EXERCISE2:
            cmd_blast = self.get_output_file(self.CMD_BLAST_LINUX, message='Where is blast executable located?')
            chosen_input_file = self.choose_input_file_menu(self.INPUT_FILES_DIRECTORY_EX2)

            if chosen_input_file == self.BACK:
                self.start()
            else:
                chosen_input_file = self.INPUT_FILES_DIRECTORY_EX2 + chosen_input_file
                output_file = self.get_output_file(self.DEFAULT_OUTPUT_EX2)
                online = self.get_menu_answer(menu=self.PROCESSING_TYPES, message="Where do you want to process?")

                if online is self.ONLINE:
                    online = True
                else:
                    online = False

                print(self.STARTED_PROCESSING)
                Exercise2.run(chosen_input_file, output_file, cmd_blast, online)
                print(self.ENDING_PROCESSING)

        elif answer == self.EXERCISE3:
            chosen_input_file = self.choose_input_file_menu(self.INPUT_FILES_DIRECTORY_EX3, message="Choose yout input MSA file")

            if chosen_input_file == self.BACK:
                self.start()
            else:
                chosen_input_file = self.INPUT_FILES_DIRECTORY_EX3 + chosen_input_file
                output_file = self.get_output_file(self.DEFAULT_OUTPUT_EX3)
                print(self.STARTED_PROCESSING)
                Exercise3.run(chosen_input_file, output_file)
                print(self.ENDING_PROCESSING)

        elif answer == self.EXERCISE4:
            chosen_input_file = self.choose_input_file_menu(self.INPUT_FILES_DIRECTORY_EX4, message="Choose yout input blast file")

            if chosen_input_file == self.BACK:
                self.start()
            else:
                chosen_input_file = self.INPUT_FILES_DIRECTORY_EX4 + chosen_input_file
                pattern = 'Mus Musculus'
                output_file = self.get_output_file(self.DEFAULT_OUTPUT_EX4 + '_' + pattern + '.out')
                print(self.STARTED_PROCESSING)
                Exercise4.run(chosen_input_file, pattern, output_file, self.EX4_RESULT_SIZE)
                print(self.ENDING_PROCESSING)

    def choose_input_file_menu(self, directory=INPUT_FILES_DIRECTORY_EX1, message="Which input file do you choose?"):
        files = self.get_files_from_directory(directory)
        files.append(self.BACK)
        return self.get_menu_answer(files, message=message)

    @staticmethod
    def get_files_from_directory(directory):
        files = [name for name in os.listdir(directory)
                 if not os.path.isdir(os.path.join(directory, name))]
        files.sort()
        return files
