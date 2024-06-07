import os
from Basic_Tools.lists_and_files import file_to_list

if __name__ == "__main__":
    cmd_file = "cmd_blast.txt"
    commands = file_to_list("cmd_blast.txt")
    for command in commands:
        os.system(command)
