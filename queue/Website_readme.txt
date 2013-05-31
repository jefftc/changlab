How to start the web2py engine in unicron

1. In unicron, go /home/xchen/web2py/
2. use screen command

-> screen 
-> python web2py.py -p 8080 -a baby2008
(ctrl-a d) optional
-> logout

3. to stop the engine, to to unicron and 
-> screen -r
ctrl-c

How to access Betsy website
1. In your local machine, set up a ssh tunnel to unicron use:
ssh -L 8080:unicron:8080 unicron

2.In your local web browser, use the 
http://localhost:8080/Betsy/default/index

Then you will access the Betsy website. 
When you first time use, you have to register an account.

