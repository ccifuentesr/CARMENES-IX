def process_bib_file(file_path):
    with open(file_path, 'r', encoding='utf-8') as file:
        lines = file.readlines()

    unique_blocks = []
    current_block = []

    for line in lines:
        if line.startswith('@'):
            # Nuevo bloque encontrado, agrega el bloque actual a la lista única
            unique_blocks.append(''.join(current_block))
            current_block = [line]
        else:
            # Agrega la línea al bloque actual
            current_block.append(line)

    # Agrega el último bloque después del bucle
    unique_blocks.append(''.join(current_block))

    # Elimina duplicados dentro de cada bloque
    unique_blocks = list(set(unique_blocks))

    with open(file_path, 'w', encoding='utf-8') as file:
        file.writelines('\n\n'.join(unique_blocks))

# Reemplaza 'biblio.bib' con la ruta correcta de tu archivo
process_bib_file('Data/biblio.bib')
