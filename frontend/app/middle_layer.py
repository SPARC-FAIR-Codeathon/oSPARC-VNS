
# the middle layer is responsible for keeping the 'task' alive
# after the user presses run, moving the SPARC api through the 
# necessary intermediate steps even if the user navigates away 
# from the front-end. The last celery task is to send an email
# to the user with their results. 

from celery import Celery

app = Celery('vnsModel',broker='redis://localhost')



@app.task
def add(x, y):
    return x + y
