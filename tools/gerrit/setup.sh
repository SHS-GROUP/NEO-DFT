#!/usr/bin/env bash

# Run this script to set up your clone of GAMESS to submit patches to Gerrit.

# Set up user name and email address
setup_git_user() {
  read -ep "Please enter your full name, e.g. 'John E. Doe': " name
  echo "Name: '$name'"
  git config user.name "$name"
  read -ep "Please enter your email address, e.g. 'john@doe.com': " email
  echo "Email address: '$email'"
  git config user.email "$email"
}

# Infinite loop until confirmation information is correct
for (( ; ; ))
do
  # Display the final user information.
  gitName=$(git config user.name)
  gitEmail=$(git config user.email)
  echo "Your commits will have the following author information:

  $gitName <$gitEmail>
"
  read -ep "Is the name and email address above correct? [Y/n] " correct
  if [ "$correct" == "n" ] || [ "$correct" == "N" ]; then
    setup_git_user
  else
    break
  fi
done

# Set up gerrit remote
gerrit_user() {
  read -ep "Enter your gerrit user (Gerrit Settings/Profile) [$USER]: " gu
  if [ "$gu" == "" ]; then
    gu=$USER
  fi
  echo -e "\nConfiguring 'gerrit' remote with user '$gu'..."
  if git config remote.gerrit.url >/dev/null; then
    # Correct the remote url
    git remote set-url gerrit ssh://$gu@www.msg.chem.iastate.edu:29418/GAMESS || \
      die "Could not set gerrit remote."
  else
    # Add a new one
    git remote add gerrit ssh://$gu@www.msg.chem.iastate.edu:29418/GAMESS || \
      die "Could not add gerrit remote."
  fi
  cat << EOF

For more information on working with Gerrit,

  http://avogadro.openmolecules.net/wiki/Working_with_Gerrit
EOF
}

# Make sure we are inside the repository.
cd "$(echo "$0"|sed 's/[^/]*$//')"

for (( ; ; ))
do
	git config remote.gerrit.url >/dev/null
	if [ $? -ne 0 ]
	then
		echo
		gerrit_user
		break
	else
		echo "The configured Gerrit remote URL is:"
		echo
		git config remote.gerrit.url
		gu=`git config remote.gerrit.url | sed -e 's/^ssh:\/\///' | sed -e 's/@www.msg.chem.iastate.edu:29418\/GAMESS//'`
		echo
		read -ep "Is the username and URL correct? [Y/n]: " correct
		if [ "$correct" == "n" ] || [ "$correct" == "N" ]; then
			gerrit_user
		else
			break
		fi
	fi
done
cat << EOF

Setting up some useful git aliases for you. This can be used by typing git and
the alias name. You can inspect all aliases in this script, or by reading
.git/config in your clone.

  prepush          - view a short form of the commits about to be pushed,
                     relative to origin/master
  gerrit-push      - push the current topic branch to Gerrit for code review.


EOF
gitroot=`git rev-parse --git-dir`
if [ ! -e ${gitroot}/hooks/commit-msg ]
then
	echo
	echo "Attempting to fetch the gerrit hook."
	echo
	echo "scp -p -P 29418 $gu@www.msg.chem.iastate.edu:hooks/commit-msg .git/hooks/commit-msg"
	scp -p -P 29418 $gu@www.msg.chem.iastate.edu:hooks/commit-msg ${gitroot}/hooks/commit-msg
	echo
fi

git config alias.prepush 'log --graph --stat origin/master..'
git_branch="\$(git symbolic-ref HEAD | sed -e 's|^refs/heads/||')"
git config alias.gerrit-push "!sh -c \"git push gerrit HEAD:refs/for/master/${git_branch}\""
