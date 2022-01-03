import sqlite3

conn = sqlite3.connect('../rm_planeten.db')
c = conn.cursor()

def read_from_db():
    c.execute('SELECT * FROM planeten')
    data = c.fetchall()
    for row in data:
        print(row[:5])

    return


read_from_db()
