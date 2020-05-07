from exercises.client_service import ClientService
from exercises.exercise1 import Exercise1
from exercises.exercise2 import Exercise2

client = ClientService()
client.start()

#por si queremos correr debugger
# Exercise1.run("archives/gb_files/NM_000492.4(homo_sapiens).gb", "output/result.fas")
#Exercise2.run("archives/prot_sequences/result.fas", "archives/blast/blast.out", online=True)


